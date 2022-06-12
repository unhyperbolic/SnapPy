from snappy.CyOpenGL import SimpleImageShaderWidget

__all__ = ['SurgeryView']

class SurgeryView(SimpleImageShaderWidget):
    fragment_shader_source = b"""
        #version 150

        out vec4 out_FragColor;
        uniform vec2 screenResolution;

        void main()
        {
            vec2 xy = (gl_FragCoord.xy - 0.5*screenResolution.xy)/screenResolution.x;
            if (xy == vec2(0,0)) {
                xy = vec2(0.001,0.001);
            }
            vec2 invz = vec2(xy.x, -xy.y)/dot(xy,xy);
            float minz = min(min(fract(invz.x), 0.1),min(fract(invz.y), 0.1));
            float maxz = min(max(1.0 - fract(invz.x), 0.1),max(1.0 - fract(invz.y), 0.1));
            vec3 col =  minz * minz * maxz * maxz * 1000.0 * vec3(1.0, 1.0, 1.0);
            out_FragColor = vec4(col, 1.0);
            if (dot(invz,invz) > 100.0) {
                out_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
            }
        }
        """
    def __init__(self,
                 container,
                 viewer,
                 cusp_index,
                 *args,
                 **kwargs):
        SimpleImageShaderWidget.__init__(
            self, container,
            *args, **kwargs)
        self.bind('<Button-1>', self.tkButton1)
        self.bind('<B1-Motion>', self.tkButtonMotion1)
        self.viewer = viewer
        self.cusp_index = cusp_index

    def get_uniform_bindings(self, width, height):
        result = {
                    'screenResolution' : ('vec2', [width, height]),
                }
        return result

    def tkButton1(self, event):
        self.process_surgery_mouse_event(event)

    def tkButtonMotion1(self, event):
        self.process_surgery_mouse_event(event)

    def process_surgery_mouse_event(self, event):
        width = self.winfo_width()
        height = self.winfo_height()
        u = (event.x - 0.5*width)/width
        v = (event.y - 0.5*height)/width
        if u*u + v*v < 0.001:
            invzu = 0
            invzv = 0
        else:
            invzu = u/(u*u + v*v)
            invzv = -v/(u*u + v*v)
        self.viewer.filling_dict['fillings'][1][self.cusp_index] = [invzu,invzv]
        self.viewer.push_fillings_to_manifold()
        self.viewer.cusp_radio_buttons[self.cusp_index].config(
                text = "Cusp %d: [%.3f,%.3f]" % (self.cusp_index, invzu, invzv))
