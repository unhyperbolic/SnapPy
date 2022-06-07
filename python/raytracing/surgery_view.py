from snappy.CyOpenGL import SimpleImageShaderWidget

__all__ = ['SurgeryView']

class SurgeryView(SimpleImageShaderWidget):
    def __init__(self,
                 container,
                 *args,
                 **kwargs):
        SimpleImageShaderWidget.__init__(
            self, container,
            *args, **kwargs)

