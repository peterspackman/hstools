"""Module for visualizing a HS with glumpy"""
import logging
import sys
from glumpy import app, gl, glm, gloo
from hstools.decompose import mean_radius, shift_to_origin
import matplotlib.pyplot as plt
import numpy as np
import sbf

LOG = logging.getLogger(__name__)

class Isosurface:
    """Wrapper class around an isosurface/mesh"""
    def __init__(self, vertices, faces, vertex_normals, vertex_colors, **kwargs):
        self.vtype = [('position', np.float32, 3),
                      ('normal', np.float32, 3),
                      ('color', np.float32, 4)]

        self.vertex_buffer = np.zeros(vertices.shape[0], self.vtype)
        self.vertex_buffer['position'] = vertices
        self.vertex_buffer['normal'] = vertex_normals
        self.vertex_buffer['color'] = vertex_colors
        self.index_buffer = np.zeros(faces.shape, dtype=np.uint32)
        self.index_buffer[:, :] = faces[:, :]
    
    @property
    def vertex_shader(self):
        return """
            uniform mat4   u_model;         // Model matrix
            uniform mat4   u_view;          // View matrix
            uniform mat4   u_projection;    // Projection matrix
            attribute vec3 position;      // Vertex position
            attribute vec3 normal;        // Vertex normal
            attribute vec4 color; // Vertex color
            varying vec3   v_normal;        // Interpolated normal (out)
            varying vec3   v_position;      // Interpolated position (out)
            varying vec4 v_color; // Vertex color

            void main()
            {
                // Assign varying variables
                v_normal   = normal;
                v_position = position;
                v_color = color;

                // Final position
                gl_Position = u_projection * u_view * u_model * vec4(position,1.0);
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
                gl_FragColor = v_color * (0.1 + 0.9*brightness * vec4(u_light_intensity, 1));
            }
        """

    @property
    def vertices(self):
        return self.vertex_buffer.view(gloo.VertexBuffer)

    @property
    def indices(self):
        return self.index_buffer.view(gloo.IndexBuffer)

    @staticmethod
    def from_sbf_file(filename, surface_property='d_norm', colormap='bwr_r'):
        surface_data = sbf.read_file(filename)
        # set center to origin
        positions = shift_to_origin(surface_data['vertices'].data.transpose())
        radius = mean_radius(positions)
        positions /= radius
        faces = surface_data['faces'].data.transpose() - 1
        vertex_normals = surface_data['vertex normals'].data.transpose()
        color_vals = surface_data[surface_property].data
        colors = values_to_colors(color_vals, colormap=plt.get_cmap(colormap))
        return Isosurface(positions, faces, vertex_normals, colors)


def values_to_colors(vals, colormap=plt.get_cmap('bwr_r')):
    maxval = np.max(vals)
    minval = np.min(vals)
    normalized = vals[:] - minval
    normalized /= (maxval - minval) 
    return colormap(normalized)


class Renderer:

    def __init__(self, render_object, **kwargs):
        self.window = app.Window(**kwargs)
        self.render_object = render_object
        self.animating = True
        self.program = gloo.Program(self.render_object.vertex_shader,
                                    self.render_object.fragment_shader)
        self.program.bind(self.render_object.vertices)
        self.program["u_light_position"] = (-2, -2, 2)
        self.program["u_light_intensity"] = (1, 1, 1)
        self.program['u_model'] = np.eye(4, dtype=np.float32)
        self.view = glm.translation(0, 0, -5) # camera position
        self.model = np.eye(4, dtype=np.float32)
        self.program['u_view'] = self.view
        self.phi, self.theta = (40, 30)

        @self.window.event
        def on_draw(dt):
            self.window.clear()
            gl.glEnable(gl.GL_DEPTH_TEST)
            self.program.draw(gl.GL_TRIANGLES, self.render_object.indices)
            # Rotate cube
            self.update_camera()

        @self.window.event
        def on_resize(width, height):
            self.program['u_projection'] = glm.perspective(45.0, width / float(height), 2.0, 100.0)

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
            if character == 'a':
                self.animating = not self.animating
            LOG.debug('Character entered (character: %s)'% character)
            LOG.debug('Character == q {}'.format(character == 'q'))

        @self.window.event
        def on_key_release(symbol, modifiers):
            LOG.debug('Key released (symbol=%s, modifiers=%s)', symbol, modifiers)

        @self.window.event
        def on_mouse_press(x, y, button):
            LOG.debug('Mouse button pressed (x=%.1f, y=%.1f, button=%d)', x, y, button)

        @self.window.event
        def on_mouse_release(x, y, button):
            LOG.debug('Mouse button released (x=%.1f, y=%.1f, button=%d)', x, y, button)

        @self.window.event
        def on_mouse_motion(x, y, dx, dy):
            LOG.debug('Mouse motion (x=%.1f, y=%.1f, dx=%.1f, dy=%.1f)',
                      x, y, dx, dy)

        @self.window.event
        def on_mouse_drag(x, y, dx, dy, button):
            LOG.debug('Mouse drag (x=%.1f, y=%.1f, dx=%.1f, dy=%.1f, button=%d)',
                      x, y, dx, dy, button)

        @self.window.event
        def on_mouse_scroll(x, y, dx, dy):
            LOG.debug('Mouse scroll (x=%.1f, y=%.1f, dx=%.1f, dy=%.1f)', x, y, dx, dy)
    
    def update_modelview(self):
        self.view = self.program['u_view'].reshape(4, 4)
        self.model = np.eye(4, dtype=np.float32)
        glm.rotate(self.model, self.theta, 0, 0, 1)
        glm.rotate(self.model, self.phi, 0, 1, 0)
        self.program['u_model'] = self.model
        self.program['u_normal'] = np.array(np.matrix(np.dot(self.view, self.model)).I.T)
    
    def update_camera(self):
        # Rotate cube
        if self.animating:
            self.theta += 0.9 # degrees
            self.phi += 1.0 # degrees
        self.update_modelview()


    def set_theta(self, value):
        self.theta = value

    def set_phi(self, value):
        self.phi = value


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
    iso = Isosurface.from_sbf_file(args.filename,
                                   surface_property=args.surface_property,
                                   colormap=args.color_map)
    renderer = Renderer(iso, width=1024, height=1024,
                        color=(0.30, 0.30, 0.35, 1.00))
    app.run(interactive=True)

if __name__ == '__main__':
    main()