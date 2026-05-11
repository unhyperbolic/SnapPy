from ..exceptions import SnapPeaFatalError

# Implementation of the SnapPea UI functions and their global variables.
cdef public void uFatalError(const_char_ptr function,
                             const_char_ptr file) except *:
    # Only raise exception the first time so that we see the first
    # uFatalError which is usually the root cause of the problem.
    if not PyErr_Occurred():
        raise SnapPeaFatalError('SnapPea crashed in function %s(), '
                                'defined in %s.c.' % (function, file))

# Global variables used for interrupt processing
cdef public Boolean gLongComputationInProgress
cdef public Boolean gLongComputationCancelled
cdef public gLongComputationTicker

# If not None, this will be called in gLongComputationContinues.
# This enables a GUI to do updates during long computations.
UI_callback = None


def SnapPea_interrupt():
    """
    The UI can call this to stop SnapPea.  Returns True if SnapPea s busy.
    If SnapPea is busy, the side effect is to set the gLongComputationCancelled
    flag.
    """
    global gLongComputationCancelled
    global gLongComputationInProgress
    if gLongComputationInProgress:
        gLongComputationCancelled = True
    return gLongComputationInProgress


cdef public void uLongComputationBegins(const_char_ptr message,
                                        Boolean is_abortable):
    global gLongComputationCancelled
    global gLongComputationInProgress
    global gLongComputationTicker
    gLongComputationCancelled = False
    gLongComputationInProgress = True
    gLongComputationTicker = time.time()

cdef public c_FuncResult uLongComputationContinues() except *:
    global gLongComputationCancelled
    global gLongComputationInProgress
    global gLongComputationTicker
    cdef now
    if gLongComputationCancelled:
        return func_cancelled
    elif UI_callback is not None:
        now = time.time()
        if now - gLongComputationTicker > 0.2:
            UI_callback()
            gLongComputationTicker = now
    return func_OK

cdef public void uLongComputationEnds() except*:
    global gLongComputationCancelled
    global gLongComputationInProgress
    gLongComputationInProgress = False
    if gLongComputationCancelled:
        gLongComputationCancelled = False
        if UI_callback is not None:
            UI_callback(interrupted=True)


show_uAcknowledge = False

cdef public void uAcknowledge(const_char_ptr message):
    if show_uAcknowledge:
        sys.stderr.write(message.decode())
        sys.stderr.write('\n')
    return

cdef public void uAbortMemoryFull():
    sys.stderr.write('Out of memory.\n')
    sys.exit(2)

cdef public int uQuery(const_char_ptr  message,
                       const_int       num_responses,
                       const_char_ptr  responses[],
                       const_int       default_response):
    #  If desired you could write this function to obtain a response
    #  from the user, but for now it is set up to return the default
    #  response, to facilitate batch computations.
    cdef char *default = <char *> responses[<int> default_response]
    sys.stderr.write('Q: %s\nA:  %s\n' % (<char *> message, default))
    return <int> default_response


