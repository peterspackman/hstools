"""Module for visualizing a HS with glumpy"""
import logging
import random
import sys
from glumpy import app, gl, glm, gloo
from glumpy.transforms import OrthographicProjection, Trackball, Position
from glumpy.graphics.text import FontManager
from glumpy.graphics.collections import GlyphCollection
from hstools.decompose import mean_radius, shift_to_origin, sht_isosurface, reconstruct_surface
import matplotlib.pyplot as plt
import numpy as np
import sbf

jab = ("`Twas brillig, and the slithy toves\n"
"  Did gyre and gimble in the wabe:\n"
"All mimsy were the borogoves,\n"
"  And the mome raths outgrabe.\n"
"\n"
"\"Beware the Jabberwock, my son!\n"
"  The jaws that bite, the claws that catch!\n"
"Beware the Jubjub bird, and shun\n"
"  The frumious Bandersnatch!\"\n"
"He took his vorpal sword in hand:\n"
"  Long time the manxome foe he sought --\n"
"So rested he by the Tumtum tree,\n"
"  And stood awhile in thought.\n"
"And, as in uffish thought he stood,\n"
"  The Jabberwock, with eyes of flame,\n"
"Came whiffling through the tulgey wood,\n"
"  And burbled as it came!\n"
"One, two! One, two! And through and through\n"
"  The vorpal blade went snicker-snack!\n"
"He left it dead, and with its head\n"
"  He went galumphing back.\n"
"\"And, has thou slain the Jabberwock?\n"
"  Come to my arms, my beamish boy!\n"
"O frabjous day! Callooh! Callay!'\n"
"  He chortled in his joy.\n"
"\n"
"`Twas brillig, and the slithy toves\n"
"  Did gyre and gimble in the wabe;\n"
"All mimsy were the borogoves,\n"
"  And the mome raths outgrabe.\n" )


LOG = logging.getLogger(__name__)

def normalize_vec3(arr):
    """ Normalize a numpy array of 3 component vectors shape=(n,3) """
    lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens
    return arr

class Isosurface:
    """Wrapper class around an isosurface/mesh"""
    def __init__(self, vertices, faces, vertex_normals, available_properties, surface_property='electric_potential'):
        self.vtype = [('position', np.float32, 3),
                      ('normal', np.float32, 3),
                      ('color', np.float32, 1)]

        self.available_properties = available_properties
        self.vertex_buffer = np.zeros(vertices.shape[0], self.vtype)
        self.vertex_buffer['position'] = vertices
        self.vertex_buffer['normal'] = vertex_normals
        self.vertex_buffer['color'] = self.available_properties[surface_property]
        self.index_buffer = np.zeros(faces.shape, dtype=np.uint32)
        self.index_buffer[:, :] = faces[:, :]
    
    @property
    def vertex_shader(self):
        return """
            uniform mat4   u_model;         // Model matrix
            uniform mat4   u_view;          // View matrix
            uniform mat4   u_projection;    // Projection matrix
            uniform mat3   u_colors;      // colors for shading
            attribute vec3 position;      // Vertex position
            attribute vec3 normal;        // Vertex normal
            attribute float color; // Vertex color
            uniform float u_min_color;
            uniform float u_max_color;
            varying vec3   v_normal;        // Interpolated normal (out)
            varying vec3   v_position;      // Interpolated position (out)
            varying vec4 v_color; // Vertex color

            vec4 colormap(float value) {
                int index = 0;
                float factor = 0;
                if (value < 0) {
                    factor = 1.0 - value / u_min_color;
                    index = 2;
                }
                else {
                    factor = 1.0 - value / u_max_color;
                }
                vec3 color = u_colors[index];
                vec3 mid = u_colors[1];
                vec3 result = color + (mid - color) * factor;
                return vec4(result[0],result[1],result[2],1);
            }

            void main()
            {
                // Assign varying variables
                v_normal   = normal;
                v_position = position;
                v_color = colormap(color);
                // Final position
                gl_Position = <transform>;
            }
        """

    @property
    def fragment_shader(self):
        return """
            uniform mat4      u_model;           // Model matrix
            uniform mat4      u_view;            // View matrix
            uniform mat4      u_normal;          // Normal matrix
            uniform mat4      u_projection;      // Projection matrix
            uniform vec3      u_light_position;  // Light position
            uniform vec3      u_light_intensity; // Light intensity

            varying vec3      v_normal;          // Interpolated normal (in)
            varying vec3      v_position;        // Interpolated position (in)
            varying vec4 v_color; // Vertex color

            void main()
            {
                // Calculate normal in world coordinates
                vec3 normal = normalize(u_normal * vec4(v_normal,1.0)).xyz;

                // Calculate the location of this fragment (pixel) in world coordinates
                vec3 position = vec3(u_view*u_model * vec4(v_position, 1));

                // Calculate the vector from this pixels surface to the light source
                vec3 surfaceToLight = u_light_position - position;

                // Calculate the cosine of the angle of incidence (brightness)
                float ref = dot(normal, surfaceToLight);
                float brightness = ref /
                                (length(surfaceToLight) * length(normal));
                brightness = max(min(brightness,1.0),0.0);

                // Calculate final color of the pixel, based on:
                // 1. The angle of incidence: brightness
                // 2. The color/intensities of the light: light.intensities

                // Final color
                gl_FragColor = v_color * (0.9 + 0.0 * brightness * vec4(u_light_intensity, 1));
            }
        """

    @property
    def vertices(self):
        return self.vertex_buffer.view(gloo.VertexBuffer)

    @property
    def indices(self):
        return self.index_buffer.view(gloo.IndexBuffer)

    def change_surface_property(self, new_property, program):
        if new_property in self.available_properties.keys():
            self.vertex_buffer['color'] = self.available_properties[new_property]

    @staticmethod
    def from_sbf_file(filename, surface_property='electric_potential',
                      offset=(0,0,0), orient=True):
        surface_data = sbf.read_file(filename)
        # set center to origin
        positions = shift_to_origin(surface_data['vertices'].data.transpose())
        radius = mean_radius(positions)
        positions /= radius * 20
        if orient:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=3)
            positions = pca.fit_transform(positions)
        offset = np.array(offset)
        positions += offset
        faces = surface_data['faces'].data.transpose() - 1
        vertex_normals = surface_data['vertex normals'].data.transpose()
        surface_properties = {}

        for dset in surface_data.datasets():
            if dset.data.shape == (positions.shape[0],):
                surface_properties[dset.name] = dset.data
        return Isosurface(positions, faces, vertex_normals, surface_properties, surface_property=surface_property)

    @staticmethod
    def from_sht_coefficients(filename, offset=(0,0,0), property='d_norm', orient=True):
        name, others, coeffs = sht_isosurface(filename, prop=property)
        mean_radius, color_min, color_scale = others
        verts, faces, colors = reconstruct_surface(coeffs, color_min=color_min,
                                                   color_scale=color_scale)

        # calculate vertex normals
        normals = np.zeros(verts.shape, dtype=verts.dtype )
        tris = verts[faces]
        n = np.cross(tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
        normalize_vec3(n)
        normals[ faces[:,0] ] += n
        normals[ faces[:,1] ] += n
        normals[ faces[:,2] ] += n
        normalize_vec3(normals)
        if orient:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=3)
            verts = pca.fit_transform(verts)

        verts /= 20
        offset = np.array(offset)
        verts += offset
        return Isosurface(verts, faces, normals, {'color': colors}, surface_property='color')


class Renderer:
    def __init__(self, render_objects, **kwargs):
        self.window = app.Window(**kwargs)
        self.render_objects = render_objects
        self.animating = True
        self.programs = []        
        for obj in self.render_objects:
            self.programs.append(gloo.Program(obj.vertex_shader,
                                 obj.fragment_shader))
            self.programs[-1].bind(obj.vertices)

        self.set_light_position([0, 0, 2])
        self.set_light_intensity([1, 1, 1])
        self.view = glm.translation(0, 0, 1) # camera position
        self.model = np.eye(4, dtype=np.float32)
        for obj, prog in zip(self.render_objects, self.programs):
            prog['u_model'] = np.eye(4, dtype=np.float32)
            prog['u_min_color'] = np.min(obj.vertex_buffer['color'])
            prog['u_max_color'] = np.max(obj.vertex_buffer['color'])
            prog['u_view'] = self.view
            prog['u_colors'] = [[0, 0, 1],
                                [1, 1, 1],
                                [1, 0, 0]]

            prog['transform'] = Trackball(Position("position"), theta=0, phi=0)
            self.window.attach(prog['transform'])
            glyphs = GlyphCollection(transform=Trackball(Position()))
            glyphs.append(jab, FontManager.get("Roboto-Regular.ttf"))
            self.window.attach(glyphs["transform"])
            self.window.attach(glyphs["viewport"])



        @self.window.event
        def on_draw(dt):
            self.window.clear()
            gl.glEnable(gl.GL_DEPTH_TEST)
            for obj, prog in zip(self.render_objects, self.programs):
                prog.draw(gl.GL_TRIANGLES, obj.indices)
            # Rotate cube
            self.update_camera()

        @self.window.timer(1 / 5.0)
        def timer(dt):
            self.window.set_title("FPS: {:3.1f}".format(self.window.fps).encode())



        @self.window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glDisable(gl.GL_BLEND)

        @self.window.event
        def on_key_press(symbol, modifiers):
            LOG.debug('Key pressed (symbol=%s, modifiers=%s)', symbol, modifiers)

        @self.window.event
        def on_character(character):
            if character == 'q':
                self.window.close()
                app.quit()
                sys.exit(0)
            if character == 'r':
                # random colormap
                for prog in self.programs:
                    prog['u_colors'] = np.random.rand(3,3) 
            if character == 'p':
                available = list(self.render_objects[0].available_properties.keys())
                self.change_surface_property(random.choice(available))

            LOG.debug('Character entered (character: %s)'% character)
            LOG.debug('Character == q {}'.format(character == 'q'))

        @self.window.event
        def on_key_release(symbol, modifiers):
            LOG.debug('Key released (symbol=%s, modifiers=%s)', symbol, modifiers)

    def update_modelview(self):
        for prog in self.programs:
            self.view = prog['u_view'].reshape(4, 4)
            self.model = np.eye(4, dtype=np.float32)
            prog['u_model'] = self.model
            prog['u_normal'] = np.array(np.matrix(np.dot(self.view, self.model)).I.T)
    
    def change_surface_property(self, new_property):

        for obj, prog in zip(self.render_objects, self.programs):
            obj.change_surface_property(new_property, prog)
            prog['u_min_color'] = np.min(obj.vertex_buffer['color'])
            prog['u_max_color'] = np.max(obj.vertex_buffer['color'])
            obj.change_surface_property(new_property, prog)
            prog.bind(obj.vertices)
    
    def update_camera(self):
        # Rotate cube
        self.update_modelview()

    def set_light_position(self, position):
        for prog in self.programs:
            prog['u_light_position'] = position

    def set_light_intensity(self, intensity):
        for prog in self.programs:
            prog["u_light_intensity"] = intensity

    def set_camera_position(self, position):
        self.view = glm.translation(*position) # camera position
        for prog in self.programs:
            prog['u_view'] = self.view
    

def nearest_square_r(n):
    return round(np.sqrt(n))

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('--color-map', '-c', default='bwr_r',
                        help='Colormap for surface')
    parser.add_argument('--surface-property', '-p', default='electric_potential',
                        help='Property to color on surface')
    parser.add_argument('--log-level', '-l', default='INFO',
                        help='Log level to print')
    parser.add_argument('--sep', '-s', default=0.2, type=float,
                        help='separation distance')
    parser.add_argument('--reconstruct', '-r', default=False, type=bool,
                        help='Reconstruct shapes via SHT')
    parser.add_argument('--orient', default=False, type=bool,
                        help='Orient the shapes via PCA')

    args = parser.parse_args()
    logging.basicConfig(level=args.log_level)
    objs = []
    n = len(args.filenames)
    if n > 1:
        col = nearest_square_r(n)
        col = col if col >= np.sqrt(n) else col + 1

        xv, yv = np.meshgrid(np.linspace(-1, 1, col),
                             np.linspace(-1, 1, col))
        for f, x, y in zip(args.filenames, xv.flatten(), yv.flatten()):
            if args.reconstruct:
                objs.append(Isosurface.from_sht_coefficients(f,
                            offset=(x*args.sep,y*args.sep, -args.sep*1.5),
                            property=args.surface_property,
                            orient=args.orient))

            objs.append(Isosurface.from_sbf_file(f,
                        surface_property=args.surface_property,
                        offset=(x*args.sep,y*args.sep, 0),
                        orient=args.orient))
    else:
        f = args.filenames[0]
        if args.reconstruct:
            objs.append(Isosurface.from_sht_coefficients(f,
                        property=args.surface_property,
                        offset=(0,0,-args.sep*1.5),
                        orient=args.orient))
        objs.append(Isosurface.from_sbf_file(f,
                    surface_property=args.surface_property,
                    orient=args.orient))

    renderer = Renderer(objs, width=1024, height=1024,
                        color=(1.0, 1.0, 1.0, 1.00))
    app.run(framerate=60)

if __name__ == '__main__':

    main()
