cdef extern from "SnapPea.h":
    ctypedef struct Real_struct:
        Real x

    ctypedef struct Complex:
        Real real
        Real imag

    Real Real_from_string(char* num_string)

    ctypedef char Boolean

    ctypedef enum c_FuncResult "FuncResult":
        func_OK = 0
        func_cancelled
        func_failed
        func_bad_input

