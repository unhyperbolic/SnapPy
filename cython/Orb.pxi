cdef extern from "parse_orb.h":
    extern c_Triangulation *orb_read_string(char *str)
