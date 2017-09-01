"""Module for visualizing a HS with glumpy"""
import logging
import random
import sys
from glumpy import app, gl, glm, gloo
from glumpy.transforms import Trackball, Position
from hstools.decompose import mean_radius, shift_to_origin
import matplotlib.pyplot as plt
import numpy as np
import sbf

LOG = logging.getLogger(__name__)

class Isosurface:
    """Wrapper class around an isosurface/mesh"""
    def __init__(self, vertices, faces, vertex_normals, available_properties, surface_property='d_norm'):
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
                float brightness = dot(normal, surfaceToLight) /
                                (length(surfaceToLight) * length(normal));
                brightness = max(min(brightness,1.0),0.0);

                // Calculate final color of the pixel, based on:
                // 1. The angle of incidence: brightness
                // 2. The color/intensities of the light: light.intensities

                // Final color
                gl_FragColor = v_color * (0.5 + 0.5 * brightness * vec4(u_light_intensity, 1));
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
    def from_sbf_file(filename, surface_property='d_norm'):
        surface_data = sbf.read_file(filename)
        # set center to origin
        positions = shift_to_origin(surface_data['vertices'].data.transpose())
        radius = mean_radius(positions)
        positions /= radius
        faces = surface_data['faces'].data.transpose() - 1
        vertex_normals = surface_data['vertex normals'].data.transpose()
        surface_properties = {}

        for dset in surface_data.datasets():
            if dset.data.shape == (positions.shape[0],):
                surface_properties[dset.name] = dset.data
        return Isosurface(positions, faces, vertex_normals, surface_properties, surface_property=surface_property)


class Renderer:
    def __init__(self, render_object, **kwargs):
        self.window = app.Window(**kwargs)
        self.render_object = render_object
        self.animating = True
        self.program = gloo.Program(self.render_object.vertex_shader,
                                    self.render_object.fragment_shader)
        self.program.bind(self.render_object.vertices)
        self.set_light_position([-2, -2, 2])
        self.set_light_intensity([1, 1, 1])
        self.program['u_model'] = np.eye(4, dtype=np.float32)
        self.program['u_min_color'] = np.min(self.render_object.vertex_buffer['color'])
        self.program['u_max_color'] = np.max(self.render_object.vertex_buffer['color'])
        self.view = glm.translation(0, 0, -5) # camera position
        self.model = np.eye(4, dtype=np.float32)
        self.program['u_view'] = self.view
        self.program['u_colors'] = [[0, 0, 1],
                                    [1, 1, 1],
                                    [1, 0, 0]]

        self.program['transform'] = Trackball(Position("position"))
        self.window.attach(self.program['transform'])

        @self.window.event
        def on_draw(dt):
            self.window.clear()
            gl.glEnable(gl.GL_DEPTH_TEST)
            self.program.draw(gl.GL_TRIANGLES, self.render_object.indices)
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
                self.program['u_colors'] = np.random.rand(3,3) 
            if character == 'p':
                available = list(self.render_object.available_properties.keys())
                self.change_surface_property(random.choice(available))

            LOG.debug('Character entered (character: %s)'% character)
            LOG.debug('Character == q {}'.format(character == 'q'))

        @self.window.event
        def on_key_release(symbol, modifiers):
            LOG.debug('Key released (symbol=%s, modifiers=%s)', symbol, modifiers)

    def update_modelview(self):
        self.view = self.program['u_view'].reshape(4, 4)
        self.model = np.eye(4, dtype=np.float32)
        self.program['u_model'] = self.model
        self.program['u_normal'] = np.array(np.matrix(np.dot(self.view, self.model)).I.T)
    
    def change_surface_property(self, new_property):
        self.render_object.change_surface_property(new_property, self.program)
        self.program['u_min_color'] = np.min(self.render_object.vertex_buffer['color'])
        self.program['u_max_color'] = np.max(self.render_object.vertex_buffer['color'])
        # rebind to the program
        self.program.bind(self.render_object.vertices)
    
    def update_camera(self):
        # Rotate cube
        self.update_modelview()

    def set_light_position(self, position):
        self.program['u_light_position'] = position

    def set_light_intensity(self, intensity):
        self.program["u_light_intensity"] = intensity

    def set_camera_position(self, position):
        self.view = glm.translation(*position) # camera position
        self.program['u_view'] = self.view
    

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', default=None)
    parser.add_argument('--color-map', '-c', default='bwr_r',
                        help='Colormap for surface')
    parser.add_argument('--surface-property', '-p', default='d_norm',
                        help='Property to color on surface')
    parser.add_argument('--log-level', '-l', default='INFO',
                        help='Log level to print')
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level)
    iso = Isosurface.from_sbf_file(args.filename, surface_property=args.surface_property)
    renderer = Renderer(iso, width=1024, height=1024,
                        color=(0.30, 0.30, 0.35, 1.00))
    app.run(framerate=60)

if __name__ == '__main__':
    main()
