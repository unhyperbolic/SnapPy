from . import number
from .pari import pari

def to_byte_str(s):
    if isinstance(s, bytes):
        return s
    return s.encode('utf-8')

def to_str(s):
    return s.decode()

# C type for a function of Real returning an object
ctypedef object (*func_real_to_obj)(Real)

# Convert Real to gen in an appropriate manner for this environment
cdef func_real_to_obj Real2gen

if hasattr(pari, '_real_coerced_to_bits_prec'):  # Cypari
    Real2gen = Real2gen_direct
else:
    Real2gen = Real2gen_string

cdef Complex2gen(Complex C):
    """
    Convert a Complex to a pari gen.
    """
    cdef real_part = Real2gen(C.real)
    cdef imag_part = Real2gen(C.imag)
    return pari.complex(real_part, imag_part)

cdef RealImag2gen(Real R, Real I):
    return pari.complex(Real2gen(R), Real2gen(I))

cdef Complex2complex(Complex C):
    """
    Convert a Complex to a python complex.
    """
    return complex( float(<double>C.real), float(<double>C.imag) )

cdef Real2float(Real R):
    """
    Convert a Real to a python float.
    """
    return float(<double>R)

cdef Complex complex2Complex(complex z):
    """
    Convert a python complex to a Complex.
    """
    cdef Complex result
    result.real = <Real>z.real
    result.imag = <Real>z.imag
    return result

cdef Real Object2Real(obj):
    cdef char* c_string
    try:
        string = obj.as_string() if isinstance(obj, Number) else str(obj)
        # Pari idiosyncratically formats small and large numbers as,
        # e.g., "1.0 E-10" (note the space before "E").
        # Remove it - otherwise it cannot be parsed.
        string = string.replace(' ', '')
        float(string)
    except:
        raise ValueError('Cannot convert %s to a Real.' % type(obj))
    string = to_byte_str(string)
    c_string = string
    return Real_from_string(c_string)

cdef Complex Object2Complex(obj):
    cdef Real real, imag
    cdef Complex result
    if hasattr(obj, 'real') and hasattr(obj, 'imag'):
        try:
            float(obj.real)
            real = Object2Real(obj.real)
        except TypeError:  # Probably Sage type
            real = Object2Real(obj.real())
        try:
            float(obj.imag)
            imag = Object2Real(obj.imag)
        except TypeError:  # Probably Sage type
            imag = Object2Real(obj.imag())
    else:
        real = Object2Real(obj)
        imag = <Real>0.0
    result.real = real
    result.imag = imag
    return result


cdef double Real2double(Real R):
    cdef double* quad = <double *>&R
    return quad[0]

cdef B2B(Boolean B):
    return B != 0
