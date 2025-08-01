from .cache import SnapPyCache
import low_index

_low_index_version = [int(n) for n in low_index.version().split('.')]

cdef class Triangulation():
    """
    A Triangulation object represents a compact 3-manifold with torus
    boundary components, given as an ideal triangulation of the
    manifold's interior.  A Dehn-filling can be specified for each
    boundary component, allowing the description of closed 3-manifolds
    and some orbifolds.  For non-orientable 3-manifolds, the boundary
    components can also be Klein bottles. Two Triangulations are equal
    ('==') if they represent combinatorially isomorphic
    triangulations.  A Triangulation does *not* have any geometric
    structure, and usually one works with the subclass Manifold which
    adds this.  Here's a quick example:

    >>> M = Triangulation('9_42')
    >>> M.num_tetrahedra()
    5
    >>> M.is_orientable()
    True

    A Triangulation can be specified in a number of ways, e.g.

    - Triangulation('9_42') : The complement of the knot 9_42 in S^3.
    - Triangulation('m125(1,2)(4,5)') : The SnapPea census manifold m125
      where the first cusp has Dehn filling (1,2) and the second cusp has
      filling (4,5).
    - Triangulation() : Opens a link editor window where can you
      specify a link complement.

    In general, the specification can be from among the below, with
    information on Dehn fillings added.

    - SnapPea cusped census manifolds: e.g. 'm123', 's123', 'v123'.

    - Link complements:

        + Rolfsen's table: e.g. '4_1', '04_1', '5^2_6', '6_4^7', 'L20935', 'l104001'.
        + Knots and links up to 14 crossings from tabulations by Hoste
          and Thistlethwaite: e.g. 'K12a456' or 'L13n579'.
        + Hoste-Thistlethwaite Knotscape table:  e.g. '11a17' or '12n345'
        + Dowker-Thistlethwaite code: e.g. 'DT:[(6,8,2,4)]', 'DT:dadbcda'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'Braid:[1,2,-3,4]'

      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - From mapping class group data using Twister:

      'Bundle(S_{1,1}, [a0, B1])' or 'Splitting(S_{1,0}, [b1, A0], [a0,B1])'

      See the help for the 'twister' module for more.

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.

    - A Regina-style isomorphism signature, such as 'dLQbcccdxwb'.

    - A string containing the contents of a SnapPea triangulation or link
      projection file.
    """
    cdef c_Triangulation* c_triangulation
    cdef _DTcode
    cdef _PDcode
    cdef _cover_info
    cdef readonly _cache
    cdef readonly LE
    cdef _link_file_full_path
    cdef hyperbolic_structure_initialized

    def __cinit__(self, spec=None, remove_finite_vertices=True):
        if UI_callback is not None:
            uLongComputationBegins('Constructing a manifold', 1)
            UI_callback()
            uLongComputationEnds()
        # Answers to potentially hard computations are cached
        self._cache = SnapPyCache()
        self._DTcode = None
        self._PDcode = None
        self._cover_info = None
        self.LE = None
        self.hyperbolic_structure_initialized = False
        self._link_file_full_path = None
        for attr in ['__snappy__', 'snapPea', '_to_string']:
            if hasattr(spec, attr):
                spec = getattr(spec, attr)()
                break
        if spec is not None and spec != 'empty':
            if not isinstance(spec, (str, bytes)):
                raise TypeError(triangulation_help %
                                self.__class__.__name__)
            self.get_triangulation(spec, remove_finite_vertices)
            if self.c_triangulation == NULL:
                raise RuntimeError('An empty triangulation was generated.')
        elif spec is None:
            self.get_from_new_plink()

        if self.c_triangulation != NULL and not self.hyperbolic_structure_initialized:
            remove_hyperbolic_structures(self.c_triangulation)

    cdef get_from_new_plink(self, file_name=None):
        if LinkEditor is None:
            raise RuntimeError('PLink was not imported.')
        self.LE = LinkEditor(no_arcs=True,
                             callback=_plink_callback,
                             cb_menu='Send to SnapPy',
                             file_name=file_name,
                             manifold=self)
        print('Starting the link editor.\n'
              'Select Tools->Send to SnapPy to load the '
              'link complement.')

    cdef get_triangulation(self, spec, remove_finite_vertices=True):
        # Step -1 Check for an entire-triangulation-file-in-a-string
        if isinstance(spec, bytes) and spec.startswith(b'pickle:'):
            return self._from_pickle(spec, remove_finite_vertices)

        if (isinstance(spec, str) and spec.startswith('% Triangulation') or
                isinstance(spec, bytes) and spec.startswith(b'% Triangulation')):
            return self._from_string(spec, remove_finite_vertices)

        # Get fillings, if any
        m = split_filling_info.match(spec)
        name = m.group(1)
        fillings = eval( '[' + m.group(2).replace(')(', '),(')+ ']', {})

        # Step 1. The easy databases
        for db in database.__all_tables__.values():
            try:
                db._one_manifold(name, self)
                break
            except KeyError:
                pass

        # Step 2. Alternate names for the Rolfsen links
        if self.c_triangulation == NULL:
            for regex in rolfsen_link_regexs:
                m = regex.match(name)
                if not m:
                    continue
                if int(m.group('components')) > 1:
                    rolfsen_name = '%d^%d_%d' % (int(m.group('crossings')),
                                                 int(m.group('components')),
                                                 int(m.group('index')))
                else:
                    rolfsen_name = '%d_%d' % (int(m.group('crossings')),
                                              int(m.group('index')))
                database.LinkExteriors._one_manifold(rolfsen_name, self)

        # Step 3. Hoste-Thistlethwaite knots
        if self.c_triangulation == NULL:
            m = is_HT_knot.match(name)
            if m:
                self.get_HT_knot(int(m.group('crossings')), m.group('alternation'),
                                 int(m.group('index')), remove_finite_vertices)

        # Step 4. Once-punctured torus bundles
        if self.c_triangulation == NULL:
            m = is_torus_bundle.match(name)
            if m:
                self.get_punctured_torus_bundle(m)

        # Step 5. (fibered) braid complements
        if self.c_triangulation == NULL:
            m = is_braid_complement.match(name)
            if m:
                word = eval(m.group(1), {})
                num_strands = max([abs(x) for x in word]) + 1
                self.set_c_triangulation(get_fibered_manifold_associated_to_braid(num_strands, word))

        # Step 6. Dowker-Thistlethwaite codes
        if self.c_triangulation == NULL:
            m = is_int_DT_exterior.match(name)
            if m:
                code = eval(m.group(1), {})
                if isinstance(code, tuple):
                    knot = spherogram.DTcodec(*code)
                elif isinstance(code, list) and isinstance(code[0], int):
                    knot = spherogram.DTcodec([tuple(code)])
                else:
                    knot = spherogram.DTcodec(code)
                klp = knot.KLPProjection()
                self.set_c_triangulation(
                    get_triangulation_from_PythonKLP(klp, remove_finite_vertices))
                self.set_name(name)
                self._set_DTcode(knot)

            m = is_alpha_DT_exterior.match(name)
            if m:
                knot = spherogram.DTcodec(m.group(1))
                klp=knot.KLPProjection()
                self.set_c_triangulation(
                    get_triangulation_from_PythonKLP(klp, remove_finite_vertices))
                self.set_name(name)
                self._set_DTcode(knot)

        # Step 7.  Bundle or splitting is given in Twister's notation
        if self.c_triangulation == NULL:
            shortened_name = name.replace(' ', '')
            mb = is_twister_bundle.match(shortened_name)
            ms = is_twister_splitting.match(shortened_name)
            if mb or ms:
                func = bundle_from_string if mb else splitting_from_string
                tri_as_string = func(shortened_name)
                self._from_string(tri_as_string, remove_finite_vertices)

        # Step 8. Regina/Burton isomorphism signatures.
        if self.c_triangulation == NULL:
            self._from_isosig(name, remove_finite_vertices)

        # Step 9. If all else fails, try to load a manifold from a file.
        if self.c_triangulation == NULL:
            self.get_from_file(name, remove_finite_vertices)

        # Set the dehn fillings
        Triangulation.dehn_fill(self, fillings)

    cdef get_HT_knot(self, crossings, alternation, index, remove_finite_vertices):
        cdef c_Triangulation* c_triangulation
        DT = [get_HT_knot_DT(crossings, alternation, index)]
        knot = spherogram.DTcodec(DT)
        c_triangulation = get_triangulation_from_PythonKLP(
            knot.KLPProjection(), remove_finite_vertices)
        name = to_byte_str('%d' % crossings + alternation + '%d' % index)
        set_triangulation_name(c_triangulation, name)
        self._set_DTcode(knot)
        self.set_c_triangulation(c_triangulation)

    cdef get_punctured_torus_bundle(self, match):
        cdef LRFactorization* gluing
        cdef size_t LRlength, i
        cdef char* LRstring

        LRpart = to_byte_str(match.group(3).upper())
        LRlength = len(LRpart)
        LRstring = LRpart
        gluing = alloc_LR_factorization(<int>LRlength)
        gluing.is_available = True
        gluing.negative_determinant = 1 if match.group(1) in ['-', 'n'] else 0
        gluing.negative_trace = 0 if match.group(2) == '+' else 1
        for i from 0 <= i < LRlength:
            gluing.LR_factors[i] = LRstring[i]
        self.set_c_triangulation(triangulate_punctured_torus_bundle(gluing))
        free_LR_factorization(gluing)

    def _get_from_link_data(self, data, remove_finite_vertices=True):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        self.c_triangulation = get_triangulation_from_PythonKLP(data, remove_finite_vertices)

    cdef get_from_file(self, name, remove_finite_vertices=True):
        try:
            locations = [os.curdir, os.environ['SNAPPEA_MANIFOLD_DIRECTORY']]
        except KeyError:
            locations = [os.curdir]
        for location in locations:
            pathname = os.path.join(location, name)
            if os.path.isfile(pathname):
                file = open(pathname, 'r')
                first_line = file.readline()[:-1]
                file.close()
                if first_line.find('% Link Projection') > -1:
                    LM = LinkManager()
                    LM._from_string(open(pathname, 'r').read())
                    klp = LM.SnapPea_KLPProjection()
                    self._link_file_full_path = os.path.abspath(pathname)
                    self._set_DTcode(spherogram.DTcodec(*LM.DT_code()))
                    self._set_PDcode(LM.PD_code())
                    self.set_c_triangulation(
                        get_triangulation_from_PythonKLP(klp, remove_finite_vertices))
                else:
                    self.set_c_triangulation(read_triangulation(to_byte_str(pathname)))

        if self.c_triangulation == NULL:
            raise IOError('The manifold file %s was not found.\n%s'%
                          (name, triangulation_help % 'Triangulation or Manifold'))
        else:
            if remove_finite_vertices:
                self._remove_finite_vertices()

    def _remove_finite_vertices(self):
        """
        Removing any finite vertices by simplification and drilling.

        >>> isosig = 'kLLLLMQkccfigghjijjlnabnwnpsii'
        >>> T = Triangulation(isosig, remove_finite_vertices=False)
        >>> T.triangulation_isosig(decorated=False) == isosig
        True
        >>> T.num_cusps(),  T._num_fake_cusps()
        (0, 1)
        >>> T._remove_finite_vertices()
        >>> T.num_cusps(),  T._num_fake_cusps()
        (1, 0)
        """
        if self.c_triangulation != NULL:
            count_cusps(self.c_triangulation)
            if get_num_fake_cusps(self.c_triangulation) > 0:
                remove_finite_vertices(self.c_triangulation)
                count_cusps(self.c_triangulation)

    def cover_info(self):
        """
        If this is a manifold or triangulation which was constructed as
        a covering space, return a dictionary describing the cover.  Otherwise
        return 0.  The dictionary keys are 'base', 'type' and 'degree'.
        """
        if self._cover_info:
            return dict(self._cover_info)

    def plink(self):
        """
        Brings up a link editor window if the manifold is stored
        as a link complement in your current session.

        >>> M = Manifold('4_1') # stored as a triangulation with a link
        >>> M.link()
        <Link: 1 comp; 4 cross>
        >>> N = Manifold('m004') # stored as a triangulation without a link
        >>> N.link() #doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        ValueError: No associated link known.

        """
        if self.LE is not None:
            self.LE.reopen()
        elif self._link_file_full_path is not None:
            self.get_from_new_plink(self._link_file_full_path)
        elif self.DT_code() is not None:
            self.get_from_new_plink()
            L = spherogram.DTcodec(*self.DT_code(flips=True)).link()
            L.view(self.LE)
        else:
            raise ValueError('No associated link known.')

    def link(self):
        """
        If the manifold is stored as a link complement in your
        current session then it returns the number of components
        and crossing of the link. To view and interact with the
        link see :meth:`spherogram.Link.view`
        and :meth:`plink <snappy.Triangulation.plink>`.
        """
        if self._PDcode is not None:
            return spherogram.Link(self._PDcode)
        elif self.DT_code() is not None:
            return spherogram.DTcodec(*self.DT_code(flips=True)).link()
        else:
            raise ValueError('No associated link known.')

    cdef set_c_triangulation(self, c_Triangulation* c_triangulation):
        self.c_triangulation = c_triangulation

    def num_cusps(self, cusp_type='all') -> int:
        """
        Return the total number of cusps.  By giving the optional argument
        'orientable' or 'nonorientable' it will only count cusps of that type.

        >>> M = Triangulation('m125')
        >>> M.num_cusps()
        2
        """
        if cusp_type == 'all':
            return get_num_cusps(self.c_triangulation)
        elif cusp_type == 'orientable':
            return get_num_or_cusps(self.c_triangulation)
        elif cusp_type == 'nonorientable':
            return get_num_nonor_cusps(self.c_triangulation)
        else:
            raise ValueError("Acceptable cusp types are "
                             "['all', 'orientable', 'nonorientable'].")

    def _num_fake_cusps(self):
        """
        Returns the number of "fake" cusps, which is typically the number
        of finite vertices.

        >>> M = Triangulation('m004(1,2)')
        >>> F = M.filled_triangulation()
        >>> F.num_cusps(),  F._num_fake_cusps()
        (0, 1)
        >>> S = Triangulation('bkaagb', remove_finite_vertices=False)
        >>> S.num_cusps(),  S._num_fake_cusps()
        (0, 2)
        """
        count_cusps(self.c_triangulation)
        return get_num_fake_cusps(self.c_triangulation)

    def orientation_cover(self):
        """
        For a non-orientable Triangulation, returns the 2-fold cover which
        is orientable.

        >>> X = Triangulation('x123')
        >>> Y = X.orientation_cover()
        >>> (X.is_orientable(), Y.is_orientable())
        (False, True)
        >>> Y
        x123~(0,0)(0,0)
        >>> Y.cover_info()['type']
        'cyclic'
        """
        if self.is_orientable():
            raise ValueError('The Triangulation is already orientable.')

        cdef c_Triangulation* cover_c_triangulation = NULL
        cdef Triangulation new_tri

        cover_c_triangulation = double_cover(self.c_triangulation)
        new_tri = self.__class__('empty')
        new_tri.set_c_triangulation(cover_c_triangulation)
        new_tri.set_name(self.name() + '~')
        new_tri._cover_info = {'base'   : self.name(),
                               'type'   : 'cyclic',
                               'degree' : 2}
        return new_tri

    def is_orientable(self) -> bool:
        """
        Return whether the underlying 3-manifold is orientable.

        >>> M = Triangulation('x124')
        >>> M.is_orientable()
        False
        """
        orientability = Orientability[get_orientability(self.c_triangulation)]
        if orientability == 'orientable':
            return True
        if orientability == 'nonorientable':
            return False
        return None

    def copy(self):
        """
        Returns a copy of the triangulation.

        >>> M = Triangulation('m125')
        >>> N = M.copy()
        """
        cdef c_Triangulation* copy_c_triangulation = NULL
        cdef Triangulation new_tri

        if self.c_triangulation is NULL:
            return self.__class__('empty')
        copy_triangulation(self.c_triangulation, &copy_c_triangulation)
        new_tri = self.__class__('empty')
        new_tri.set_c_triangulation(copy_c_triangulation)
        return new_tri

    def randomize(self, blowup_multiple=4, passes_at_fours=6):
        """
        Perform random Pachner moves on the underlying triangulation,
        including some initial 3 -> 2 moves that increase the number of
        tetrahedra by blowup_multiple.

        >>> M = Triangulation('Braid:[1,2,-3,-3,1,2]')
        >>> M.randomize()
        """
        if self.c_triangulation is NULL:
            return
        randomize_triangulation_with_options(self.c_triangulation, passes_at_fours, blowup_multiple)
        self._cache.clear(message='randomize')

    def simplify(self, passes_at_fours=6):
        """
        Try to simplify the triangulation by doing Pachner moves.

        >>> M = Triangulation('12n123')
        >>> M.simplify()

        It does four kinds of moves that reduce the number of
        tetrahedra:

        * 3 -> 2 and 2 -> 0 Pacher moves, which eliminate one or two
          tetrahedra respectively.

        * On suitable valence-1 edges, does a 2 -> 3 and then 2 -> 0 move,
          which removes a tetrahedron and creates a new valence-1 edge.

        * When a 2-simplex has two edges of valence-4 giving rise to the
          suspension of a pentagon, replace these 6 tetrahedra with a
          single edge of valence 5.

        It also does random 4 -> 4 moves in hopes of setting up a
        simplfication.  The argument passes_at_fours is the number of
        times it goes through the valence-4 edges without progress
        before giving up.

        """
        if self.c_triangulation is NULL:
            return
        basic_simplification_with_options(self.c_triangulation, passes_at_fours)
        self._cache.clear(message='simplify')

    def _two_to_three(self, tet_num, face_index):
        """
        Performs a 2-3 move which removes a given face.

        The face is specified by giving the index of one of the
        tetrahedra which contains the face, as well as the index of that
        face within the tetrahedron (i.e. the index of the vertex
        opposite the face).

        If the face is not adjacent to two distinct tetrahedra, this
        function does nothing and returns a non-zero value.
        """

        cdef c_FuncResult result
        cdef c_Tetrahedron* tet

        n = valid_index(
            tet_num, self.num_tetrahedra(),
            "The specified tetrahedron (%s) does not exist.")

        f = valid_index(
            face_index, 4,
            "The specified face index (%s) is invalid.")

        tet = self.c_triangulation.tet_list_begin.next
        for i in range(n):
            tet = tet.next
        result = two_to_three(tet, f, &self.c_triangulation.num_tetrahedra)

        if result == func_OK:
            self._cache.clear(message = '2-3 move')
            self._polish_hyperbolic_structures()

        return result

    def _three_to_two(self, tet_num, edge_index):
        r"""
        Perform a 3-2 move which removes a given 3-valent edge.

        The edge is specified by giving the index of one of the tetrahedra which
        contains the edge, as well as the index of that edge within the
        tetrahedron (see below).

        If specified edge is not 3-valent or the three adjacent tetrahedra are
        not distinct, the function does nothing and returns a non-zero value.

                 1
                /|\
               / | \
              2  5  1
             /   |   \
            2--0-|----3
             \   |   /
              4  |  3
               \ | /
                \|/
                 0
        """
        cdef c_FuncResult result
        cdef c_Tetrahedron* tet
        cdef EdgeClass* where_to_resume

        _ = valid_index(
            tet_num, self.num_tetrahedra(),
            "The specified tetrahedron (%s) does not exist.")

        e = valid_index(
            edge_index, 6,
            "The specified edge index (%s) is invalid.")

        tet = self.c_triangulation.tet_list_begin.next
        for i in range(tet_num):
            tet = tet.next

        if tet.edge_class[e].order != 3:
            return func_failed

        result = three_to_two(tet.edge_class[e], &where_to_resume,
                              &self.c_triangulation.num_tetrahedra)

        if result == func_OK:
            self._cache.clear(message = '3-2 move')
            self._polish_hyperbolic_structures()

        return result

    def with_hyperbolic_structure(self):
        """
        Add a (possibly degenerate) hyperbolic structure, turning the
        :class:`Triangulation` into a :class:`Manifold`.

        >>> M = Triangulation('m004')
        >>> N = M.with_hyperbolic_structure()
        >>> N.volume() # doctest: +NUMERIC6
        2.02988321
        """
        return Manifold_from_Triangulation(self)

    def _polish_hyperbolic_structures(self):
        """
        This function does nothing for Triangulation. It is overloaded in
        Manifold.
        """
        pass

    def _empty_save(self):
        """
        Called by save when no file name is specified, so that in
        theory it can be overloaded by the UI.
        """
        if asksaveasfile:
            savefile = asksaveasfile(
                mode='w', title='Save Triangulation', defaultextension='.tri',
                filetypes = [
                    ('Triangulation and text files', '*.tri *.txt', 'TEXT'),
                    ('All text files', '', 'TEXT'),
                    ('All files', '')])
            if savefile:
                filename = savefile.name
                savefile.close()
                self.save(filename)
                return

        raise ValueError('Please specify a file name.')

    def save(self, file_name=None):
        """
        Save the triangulation as a SnapPea triangulation file.

        >>> M = Triangulation('m004')
        >>> M.save('fig-eight.tri')     #doctest: +SKIP

        To retrieve a SnapPea triangulation from the saved file
        you can do the following. The first command creates a cusped
        manifold M. The second one creates the filled manifold M1
        with Dehn coefficients (2,3).

        >>> M = Manifold('fig-eight.tri')   #doctest: +SKIP
        >>> M1 = Manifold('fig-eight.tri(2,3)')   #doctest: +SKIP

        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if file_name is None:
            self._empty_save()
        else:
            b_file_name = to_byte_str(file_name)
            write_triangulation(self.c_triangulation, b_file_name)

    def _to_string(self):
        """
        Return a string containing the contents of a SnapPea triangulation
        file.
        >>> M = Manifold('7_4')
        >>> seed = M._to_string()
        >>> N = Manifold(seed)
        >>> N == M
        True
        """
        cdef char *c_string
        cdef result
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        else:
            try:
                c_string = string_triangulation(self.c_triangulation)
                result = c_string
            finally:
                free(c_string)
            return to_str(result)

    def _from_string(self, string, remove_finite_vertices=True):
        """
        WARNING: Users should not use this function directly.  To
        create a Triangulation or Manifold or ManifoldHP from a string
        containing the contents of a triangulation file, simply do:

        >>> M = Manifold('7_4')
        >>> seed = M._to_string()
        >>> N = Manifold(seed)
        >>> N == M
        True

        Fill an empty Triangulation from a string generated by
        _to_string.
        """
        cdef c_Triangulation* c_triangulation = NULL
        if self.c_triangulation is not NULL:
            raise ValueError('The Triangulation must be empty.')
        b_string = to_byte_str(string)
        c_triangulation = read_triangulation_from_string(b_string)
        self.set_c_triangulation(c_triangulation)
        if remove_finite_vertices:
            self._remove_finite_vertices()

    def _to_bytes(self):
        """
        Return a reasonably short byte sequence which encodes the
        combinatorics of this triangulation.
        >>> M = Manifold('m125')
        >>> seed = M._to_bytes()
        >>> len(seed)
        12
        >>> N = Manifold('empty')
        >>> N._from_bytes(seed)
        >>> N == M
        True
        """
        cdef TerseTriangulation* c_terse
        cdef int n, byte
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if False in [ c.is_complete for c in self.cusp_info()]:
            raise ValueError('Byte encoding requires complete cusps.')
        c_terse = tri_to_terse(self.c_triangulation)
        byteseq = bytearray([c_terse.num_tetrahedra])
        byte, bit = 0, 0
        for n from 0 <= n < 2*c_terse.num_tetrahedra:
            if c_terse.glues_to_old_tet[n]:
                byte |= 1 << bit
            bit += 1
            if bit % 8 == 0:
                byteseq.append(byte)
                byte, bit = 0, 0
        if bit:
            byteseq.append(byte)
        for n from 0 <= n < 1 + c_terse.num_tetrahedra:
            byteseq.append(c_terse.which_old_tet[n])
        for n from 0 <= n < 1 + c_terse.num_tetrahedra:
            byteseq.append(c_terse.which_gluing[n])
        free_terse_triangulation(c_terse)
        return bytearray_to_bytes(byteseq)

    def _from_bytes(self, bytestring, *args):
        """
        Fill an empty triangulation from a byte sequence generated by
        _to_bytes.
        >>> M = Manifold('m125')
        >>> seed = M._to_bytes()
        >>> len(seed)
        12
        >>> N = Manifold('empty')
        >>> N._from_bytes(seed)
        >>> N == M
        True
        """
        cdef c_Triangulation* c_triangulation = NULL
        if self.c_triangulation is not NULL:
            raise ValueError('The Triangulation must be empty.')
        c_triangulation = triangulation_from_bytes(bytestring)
        self.set_c_triangulation(c_triangulation)

    def _from_pickle(self, bytestring, remove_finite_vertices=True):
        """
        """
        cdef c_Triangulation* c_triangulation = NULL
        if self.c_triangulation is not NULL:
            raise ValueError('The Triangulation must be empty.')
        c_triangulation = unpickle_triangulation(bytestring)
        self.set_c_triangulation(c_triangulation)
        if remove_finite_vertices:
            self._remove_finite_vertices()

    def _from_tetrahedra_gluing_data(self, tetrahedra_data,
                                     num_or_cusps=0,
                                     num_nonor_cusps=0,
                                     cusp_indices=None,
                                     peripheral_curves=None,
                                     remove_finite_vertices=True):
        """
        >>> data = [([1, 1, 2, 2], [[0, 1, 3, 2], [2, 3, 1, 0], [0, 1, 3, 2], [3, 0, 1, 2]]),
        ...         ([0, 2, 2, 0], [[0, 1, 3, 2], [1, 0, 2, 3], [2, 3, 1, 0], [3, 2, 0, 1]]),
        ...         ([1, 1, 0, 0], [[1, 0, 2, 3], [3, 2, 0, 1], [1, 2, 3, 0], [0, 1, 3, 2]])]
        >>> M = Triangulation('empty')
        >>> M._from_tetrahedra_gluing_data(data)
        >>> N = Triangulation('m007')
        >>> M == N
        True

        Check that cusp indices are transmitted. This manfiold has two
        cusps and trivial symmetry group and we interchange the two
        cusps.

        >>> M = Manifold('o9_44133')
        >>> tet_data = M._get_tetrahedra_gluing_data()
        >>> cusp_data = M._get_cusp_indices_and_peripheral_curve_data()[0]
        >>> cusp_data = [[1 - c for c in tet_cusps] for tet_cusps in cusp_data]
        >>> N = Manifold('empty')
        >>> N._from_tetrahedra_gluing_data(tet_data, num_or_cusps=2, cusp_indices=cusp_data)
        >>> M.is_isometric_to(N, True)[0].cusp_images()
        [1, 0]

        A nonorientable example.

        >>> Y = Triangulation('y885')
        >>> Y.set_peripheral_curves([([1, 0], [0, 1]), ([-1, 0], [0, -1]), ([1, 7], [1, 8])])
        >>> tets = Y._get_tetrahedra_gluing_data()
        >>> cusps, curves = Y._get_cusp_indices_and_peripheral_curve_data()
        >>> max(sum(curves, []))
        9
        >>> X = Triangulation('empty')
        >>> X._from_tetrahedra_gluing_data(tets, 1, 2, cusps, curves)
        >>> X._get_cusp_indices_and_peripheral_curve_data() == (cusps, curves)
        True
        """
        cdef c_Triangulation* c_triangulation = NULL
        if self.c_triangulation is not NULL:
            raise ValueError('The Triangulation must be empty.')
        c_triangulation = listlike_to_triangulation(tetrahedra_data,
                                                    num_or_cusps,
                                                    num_nonor_cusps,
                                                    cusp_indices,
                                                    peripheral_curves)
        self.set_c_triangulation(c_triangulation)
        if remove_finite_vertices:
            self._remove_finite_vertices()

    def _from_isosig(self, isosig, remove_finite_vertices=True):
        """
        WARNING: Users should not use this function directly.  To
        create a Triangulation or Manifold or ManifoldHP from an isosig,
        simply do:

        >>> M = Manifold('empty')
        >>> M._from_isosig('jLAwLMQbcbdghgiiixxnxhhxkwv')
        >>> N = Manifold('o9_10894')
        >>> M == N
        True

        Fill an empty Triangulation from an isosig generated by
        triangulation_isosig.
        """

        if self.c_triangulation is not NULL:
            raise ValueError('The Triangulation must be empty.')

        match = is_decorated_isosig.match(isosig)
        if match:
            isosig, decoration = match.groups()
        elif is_isosig.match(isosig):
            decoration = None
        else:
            # Did not match regular expression, bail
            return

        self.set_c_triangulation(
            triangulation_from_isomorphism_signature(isosig.encode('ascii')))

        # Not a valid isomorphism signature, bail
        if self.c_triangulation == NULL:
            return

        if decoration:
            decorated_isosig.set_peripheral_from_decoration(self, decoration)

        if remove_finite_vertices:
            self._remove_finite_vertices()



    def __reduce__(self):
        """
        Used by pickle.dumps.

        >>> from pickle import loads, dumps
        >>> from spherogram import random_link
        >>> M = Manifold('m125')
        >>> M == loads(dumps(M))
        True
        >>> M.dehn_fill((2,3),0)
        >>> len(dumps(M)) < len(M._to_string())
        True
        >>> N = loads(dumps(M))
        >>> N
        m125(2,3)(0,0)
        >>> M.is_isometric_to(N)
        True
        >>> M = Manifold('m024')
        >>> M.is_orientable()
        False
        >>> N = loads(dumps(M))
        >>> N == M
        True
        >>> N.volume() == M.volume()
        True
        >>> while True:
        ...    M = random_link(500).exterior()
        ...    if M.num_tetrahedra() > 256:
        ...       break
        ...
        >>> M == loads(dumps(M))
        True
        """
        return (self.__class__, (self.pickle(),))

    def pickle(self):
        try:
            return pickle_triangulation(self.c_triangulation)
        except ValueError:
            # When M has > 255 cusps or does not have integer filling
            # coefficients and curve weights in the interval [-128, 127] we use
            # _to_string.
            return self._to_string()

    def _reindex_cusps(self, permutation):
        """
        >>> import itertools
        >>> M = Manifold('L14n62484')
        >>> ['%.5f' % z.real() for z in M.cusp_info('shape')]
        ['0.10384', '1.17955', '-1.89674', '0.86324']
        >>> perms = itertools.permutations(range(4))
        >>> for perm in perms:
        ...     shapes = M.cusp_info('shape')
        ...     M._reindex_cusps(perm)
        ...     new_shapes = M.cusp_info('shape')
        ...     for i in range(4):
        ...         assert new_shapes[perm[i]] == shapes[i]
        """
        cdef int* indices
        cdef int n, num = self.num_cusps()
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if ( len(permutation) != num or set(permutation) != set(range(num)) ):
            raise ValueError('Not a valid permutation')
        indices = <int *>malloc(num*sizeof(int))
        for n in range(num):
            indices[n] = permutation[n]
        reindex_cusps(self.c_triangulation, indices)
        free(indices)

    def _get_cusp_indices_and_peripheral_curve_data(self):
        cdef int i, j, k, v, f
        cdef TriangulationData* data
        triangulation_to_data(self.c_triangulation, &data)

        result_cusp_indices = []
        for i in range(self.num_tetrahedra()):
            row = [data.tetrahedron_data[i].cusp_index[v]
                   for v in range(4)]
            result_cusp_indices.append(row)

        result_curves = []
        for i in range(self.num_tetrahedra()):
            for j in range(2):         # meridian, longitude
                for k in range(2):     # righthanded, lefthanded
                    data_ijk = data.tetrahedron_data[i].curve[j][k]
                    row = [data_ijk[v][f] for v in range(4) for f in range(4)]
                    result_curves.append(row)

        free_triangulation_data(data)
        return (result_cusp_indices, result_curves)

    def _get_tetrahedra_gluing_data(self):
        cdef int i, j, k
        cdef TriangulationData* data
        triangulation_to_data(self.c_triangulation, &data)
        result = []
        for i from 0 <= i < data.num_tetrahedra:
            neighbors = [data.tetrahedron_data[i].neighbor_index[j]
                         for j in range(4)]
            perms = [[data.tetrahedron_data[i].gluing[j][k]
                      for k in range(4)] for j in range(4)]

            result.append((neighbors, perms))
        free_triangulation_data(data)
        return result

    def _isomorphisms_to(self, Triangulation other not None):
        """
        Returns a complete list of combinatorial isomorphisms between
        the two triangulations:

        >>> M = Manifold('5^2_1')
        >>> N = Manifold('5^2_1')
        >>> N.set_peripheral_curves([[[2,3],[-1,-1]],[[1,1],[0,1]]])
        >>> isoms = M.isomorphisms_to(N)
        >>> isoms[6]
        0 -> 1  1 -> 0
        [ 1 0]  [-1 1]
        [-1 1]  [-3 2]
        Does not extend to link

        Each transformation between cusps is given by a matrix which
        acts on the left.  That is, the two *columns* of the matrix
        give the image of the meridian and longitude respectively.  In
        the above example, the meridian of cusp 0 is sent to the
        meridian of cusp 1.
        """
        cdef IsometryList *isometries = NULL

        if self.c_triangulation is NULL or other.c_triangulation is NULL:
            raise ValueError('Manifolds must be non-empty.')

        compute_cusped_isomorphisms(self.c_triangulation,
                                    other.c_triangulation,
                                    &isometries,
                                    NULL)

        if isometry_list_size(isometries) == 0:
            result = []
        else:
            result = IsometryListToIsometries(isometries)
        free_isometry_list(isometries)
        return result

    def isomorphisms_to(self, Triangulation other not None):
        """
        This is only here to make the tests that use SnapPy.Triangulation work.
        """
        return self._isomorphisms_to(other)

    def __dealloc__(self):
        if self.c_triangulation is not NULL:
            free_triangulation(self.c_triangulation)

    def __richcmp__(Triangulation self, Triangulation other, op):
        """
        Two triangulations are equal if they are combinatorially
        isomorphic.  Currently we don't handle the case where there
        are non-trivial Dehn fillings.

        >>> M = Triangulation('m004')
        >>> N = M.copy()
        >>> N == M
        True
        >>> M.dehn_fill( (5,3), 0)
        >>> N == M
        ValueError: Can't compare triangulations of manifolds with Dehn fillings.
        """
        if op != 2:
            return NotImplemented
        if type(self) != type(other):
            return False
        if False in ( self.cusp_info('is_complete') +
                      other.cusp_info('is_complete') ):
            raise ValueError("Can't compare triangulations of manifolds "
                             "with Dehn fillings.")
        return bool(same_triangulation(self.c_triangulation,
                                       other.c_triangulation))

    def __repr__(self):
        if self.c_triangulation is NULL:
            return 'Empty Triangulation'
        else:
            repr = self.name()
            for i in range(self.num_cusps()):
                info = self.cusp_info(i)
                if info.is_complete:
                    repr += '(0,0)'
                else:
                    repr += '(%g,%g)'% info['filling']
            return repr

    def name(self) -> str:
        """
        Return the name of the triangulation.

        >>> M = Triangulation('4_1')
        >>> M.name()
        '4_1'
        """
        if self.c_triangulation is NULL:
            return
        return to_str(get_triangulation_name(self.c_triangulation))

    def set_name(self, new_name : str) -> None:
        """
        Give the triangulation a new name.

        >>> M = Triangulation('4_1')
        >>> M.set_name('figure-eight-comp')
        >>> M
        figure-eight-comp(0,0)
        """
        b_new_name = to_byte_str(new_name)
        cdef char* c_new_name = b_new_name
        if self.c_triangulation is NULL:
            raise ValueError('The empty triangulation has no name.')
        set_triangulation_name(self.c_triangulation, c_new_name)

    def DT_code(self, alpha=False, flips=False):
        """
        Return the Dowker-Thistlethwaite code of this link complement,
        if it is a link complement. The DT code is intended to be an
        immutable attribute, for use with knot and link exteriors
        only, which is set only when the manifold was created.

        Here is the Whitehead link:

        >>> M = Manifold('L5a1')
        >>> M.DT_code()
        [(6, 8), (2, 10, 4)]
        >>> M.DT_code(alpha=True)
        'ebbccdaeb'
        >>> M.DT_code(alpha=True, flips=True)
        'ebbccdaeb.01110'
        >>> M.DT_code(flips=True)
        ([(6, 8), (2, 10, 4)], [0, 1, 1, 1, 0])
        """
        codec = self._DTcode
        if codec is None:
            return None
        elif not isinstance(codec, spherogram.DTcodec):
            self._DTcode = codec = spherogram.DTcodec(codec)
        if alpha:
            return codec.encode(header=False, flips=flips)
        else:
            if flips:
                return codec.code, [int(x) for x in codec.flips]
            else:
                return codec.code

    def _set_DTcode(self, code, lazy=True):
        if not lazy and not isinstance(code, spherogram.DTcodec):
            code = spherogram.DTcodec(code)
        self._DTcode = code

    def _set_PDcode(self, code):
        self._PDcode = code

    def num_tetrahedra(self) -> int:
        """
        Return the number of tetrahedra in the triangulation.

        >>> M = Triangulation('m004')
        >>> M.num_tetrahedra()
        2
        """
        if self.c_triangulation is NULL:
            return 0
        return get_num_tetrahedra(self.c_triangulation)

    def dehn_fill(self, filling_data, which_cusp=None) -> None:
        """
        Set the Dehn filling coefficients of the cusps.  This can be
        specified in the following ways, where the cusps are numbered
        by 0,1,...,(num_cusps - 1).

        - Fill cusp 2:

          >>> M = Triangulation('8^4_1')
          >>> M.dehn_fill((2,3), 2)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(0,0)

        - Fill the last cusp:

          >>> M.dehn_fill((1,5), -1)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(1,5)

        - Fill the first two cusps:

          >>> M.dehn_fill( [ (3,0), (1, -4) ])
          >>> M
          8^4_1(3,0)(1,-4)(2,3)(1,5)

        - When there is only one cusp, there's a shortcut

          >>> N = Triangulation('m004')
          >>> N.dehn_fill( (-3,4) )
          >>> N
          m004(-3,4)

        Does not return a new :class:`Triangulation`.
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        num_cusps = self.num_cusps()

        if which_cusp is not None:
            which_cusp = valid_index(
                which_cusp, num_cusps,
                'The specified cusp (%s) does not exist.')

            meridian, longitude = filling_data
            complete = ( meridian == 0 and longitude == 0)
            set_cusp_info(self.c_triangulation,
                          which_cusp, complete,
                          Object2Real(meridian),
                          Object2Real(longitude))
            self._cache.clear(message='dehn_fill')
        else:
            if num_cusps > 1 and len(filling_data) == 2:
                if ( not hasattr(filling_data, '__getitem__')
                     or not hasattr(filling_data[0], '__getitem__') ):
                    raise IndexError('If there is more than one cusp '
                                     'you must specify which one you\n'
                                     'are filling, e.g. M.dehn_fill((2,3),1)')
            if num_cusps == 1 and len(filling_data) == 2:
                Triangulation.dehn_fill(self, filling_data, 0)
                return
            if len(filling_data) > num_cusps:
                raise IndexError('You provided filling data for too '
                                 'many cusps.  There are only %s.'%
                                 num_cusps)
            for i, fill in enumerate(filling_data):
                Triangulation.dehn_fill(self, fill, i)

    # When doctesting, the M,L coefficients acquire an accuracy of 8.
    # So we have to include the zeros in the doctest string.
    def cusp_info(self, data_spec=None):
        """
        Returns an info object containing information about the given
        cusp.   Usage:

        >>> M = Triangulation('v3227(0,0)(1,2)(3,2)')
        >>> M.cusp_info(1)
        Cusp 1 : torus cusp with Dehn filling coefficients (M, L) = (1.0, 2.0)
        >>> c = M.cusp_info(1)
        >>> c.is_complete
        False
        >>> sorted(c.keys())
        ['filling', 'index', 'is_complete', 'topology']

        You can get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : torus cusp, not filled,
         Cusp 1 : torus cusp with Dehn filling coefficients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coefficients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('is_complete')
        [True, False, False]
        """
        cdef int cusp_index
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef Real m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_accuracy, current_shape_accuracy,
        cdef Complex initial_modulus, current_modulus

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if data_spec is None:
            return ListOnePerLine([self.cusp_info(i)
                                   for i in range(self.num_cusps())])
        if isinstance(data_spec, str):
            return [c[data_spec] for c in self.cusp_info()]
        cusp_index = valid_index(
            data_spec, self.num_cusps(),
            'The specified cusp (%s) does not exist.')

        get_cusp_info(self.c_triangulation, cusp_index,
                      &topology, &is_complete, &m, &l,
                      &initial_shape, &current_shape,
                      &initial_shape_accuracy, &current_shape_accuracy,
                      &initial_modulus, &current_modulus)
        info = {
           'index' : cusp_index,
           'topology' : CuspTopology[topology],
           'is_complete' : B2B(is_complete),
           'filling' : (Real2float(m), Real2float(l))
           }

        return CuspInfo(**info)

    def reverse_orientation(self) -> None:
        """
        Reverses the orientation of the Triangulation, presuming that
        it is orientable.

        >>> M = Manifold('m015')
        >>> cs = M.chern_simons()
        >>> M.reverse_orientation()
        >>> abs(cs + M.chern_simons()) # doctest: +NUMERIC9
        0.0
        """
        if not self.is_orientable():
            raise ValueError("The Manifold is not orientable, so its "
                             "orientation can't be reversed.")
        reorient(self.c_triangulation)
        self._cache.clear(message='reverse_orientation')

    def filled_triangulation(
            self, cusps_to_fill='all') -> 'snappy.Triangulation':
        """
        Return a new manifold where the specified cusps have been
        permanently filled in.  Examples:

        Filling all the cusps:

        >>> M = Triangulation('m125(1,2)(3,4)')
        >>> N = M.filled_triangulation()
        >>> N.num_cusps()
        0

        Filling cusps 0 and 2 :

        >>> M = Triangulation('v3227(1,2)(3,4)(5,6)')
        >>> M.filled_triangulation([0,2])
        v3227_filled(3,4)
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        n = self.num_cusps()
        if cusps_to_fill == 'all':
            cusps_to_fill = [c for c in range(n) if cusp_is_fillable(self.c_triangulation, c)]

        if False in [(c in range(n)) for c in cusps_to_fill]:
            raise IndexError('The specified indices to be filled are beyond '
                             'the actual number of cusps.')
        if 0 in [cusp_is_fillable(self.c_triangulation, c)
                 for c in cusps_to_fill]:
            raise IndexError('To permanently fill a cusp, the Dehn '
                             'filling coefficients must be relatively '
                             'prime integers.')

        cdef c_Triangulation* c_filled_tri = NULL
        cdef Triangulation filled_tri
        cdef Boolean *fill_cusp_spec = NULL

        fill_cusp_spec = <Boolean*>malloc(n*sizeof(Boolean))
        for i in range(n):
            fill_cusp_spec[i] = True if i in cusps_to_fill else False
        fill_all = all([i in cusps_to_fill for i in range(n)])
        c_filled_tri = fill_cusps(self.c_triangulation,
                                  fill_cusp_spec, '', fill_all)
        free(fill_cusp_spec)
        filled_tri = self.__class__('empty')
        filled_tri.set_c_triangulation(c_filled_tri)
        filled_tri.set_name(self.name() + '_filled')
        return filled_tri

    def _unsimplified_filled_triangulation(self, method='fold'):
        """
        For a Triangulation that describes a closed manifold, returns
        the unsimplified finite triangulation that the kernel builds.

        >>> M = Triangulation('m004(1, 2)')
        >>> F = M._unsimplified_filled_triangulation()
        >>> F.num_tetrahedra(), F._num_fake_cusps()
        (58, 7)

        The default is to use the kernel's original fold method as
        the final step to close off the cusp. Specifying
        method='layered' uses a 1-tetrahedron solid torus instead.
        For method='layered_and_marked', the core solid tori are marked
        internally so that any calls to simplify/randomize will leave
        these in place.

        >>> M = Triangulation('m004(1, 2)')
        >>> F = M._unsimplified_filled_triangulation(method='layered')
        >>> F.num_tetrahedra(), F._num_fake_cusps()
        (60, 7)
        >>> F.simplify(); F.num_tetrahedra() < 15
        True
        >>> M = Triangulation('K8n1(1,0)')
        >>> F = M._unsimplified_filled_triangulation(method='layered_and_marked')
        >>> F.num_tetrahedra() > 2
        True

        You can determine which tets are the cores, and remove the
        marks to allow unrestricted simplification.

        >>> any(F._marked_tetrahedra(clear_marks=True))
        True
        >>> F.simplify()
        >>> F.num_tetrahedra() <= 2
        True
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        n = self.num_cusps()
        if not all(cusp_is_fillable(self.c_triangulation, c) for c in range(n)):
            raise IndexError('All cusps must be fillable.')
        cdef c_Triangulation* c_new_tri = NULL
        cdef Triangulation filled_tri
        cdef Boolean *fill_cusp_spec = NULL
        cdef Boolean fill_by_fold
        if method == 'fold':
            fill_by_fold = True
        elif method == 'layered':
            fill_by_fold = False
            mark_solid_tori = False
        elif method == 'layered_and_marked':
            fill_by_fold = False
            mark_solid_tori = True
        else:
            raise ValueError("The method must be 'fold' or 'layered' or 'layered_and_marked'")

        c_new_tri = subdivide(self.c_triangulation, to_byte_str(self.name() + '_filled'))
        fill_cusp_spec = <Boolean*>malloc(n*sizeof(Boolean))
        for i in range(n):
            fill_cusp_spec[i] = True

        close_cusps(c_new_tri, fill_cusp_spec, fill_by_fold, mark_solid_tori)
        number_the_tetrahedra(c_new_tri)
        number_the_edge_classes(c_new_tri)
        create_fake_cusps(c_new_tri)
        count_cusps(c_new_tri)
        filled_tri = _triangulation_class('empty')
        filled_tri.set_c_triangulation(c_new_tri)
        free(fill_cusp_spec)
        return filled_tri

    def _marked_tetrahedra(self, clear_marks=False):
        cdef int *marked = NULL
        n = self.num_tetrahedra()
        marked = <int*>malloc(n*sizeof(int))
        unchangeable_tetrahedra(self.c_triangulation, marked)
        ans = [marked[i] for i in range(n)]
        free(marked)
        if clear_marks:
            all_tetrahedra_changeable(self.c_triangulation)
        return ans

    def edge_valences(self):
        """
        Returns a dictionary whose keys are the valences of the edges
        in the triangulation, and the value associated to a key is the
        number of edges of that valence.

        >>> M = Triangulation('v3227')
        >>> M.edge_valences()     # doctest: +SKIP
        {10: 1, 4: 1, 5: 2, 6: 3}
        """
        cdef int c, v = 1
        ans = {}
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        while get_num_edge_classes(self.c_triangulation, v, 1) > 0:
            c = get_num_edge_classes(self.c_triangulation, v, 0)
            if c > 0:
                ans[v] = c
            v += 1
        return ans

    def gluing_equations_pgl(self, N =2, equation_type='all'):

        """
        Returns a NeumannZagierTypeEquations object that contains a matrix
        encoding the gluing equations for boundary-parabolic PGL(N,C)
        representations together with explanations of the meaning
        of the rows and the columns of the matrix.

        This method generalizes gluing_equations() to PGL(N,C)-representations
        as described in
        Stavros Garoufalidis, Matthias Goerner, Christian K. Zickert:
        "Gluing Equations for PGL(n,C)-Representations of 3-Manifolds"
        (http://arxiv.org/abs/1207.6711).

        The result of the :meth:`gluing_equations` can be obtained from
        the general method by:

        >>> M = Triangulation('m004')
        >>> M.gluing_equations_pgl().matrix
        [ 2  1  0  1  0  2]
        [ 0  1  2  1  2  0]
        [ 1  0  0  0 -1  0]
        [ 0  0  0  0 -2  2]

        But besides the matrix, the method also returns explanations of
        the columns and rows:

        >>> M = Triangulation("m004")
        >>> M.gluing_equations_pgl()
        NeumannZagierTypeEquations(
          [ 2  1  0  1  0  2]
          [ 0  1  2  1  2  0]
          [ 1  0  0  0 -1  0]
          [ 0  0  0  0 -2  2],
          explain_columns = ['z_0000_0', 'zp_0000_0', 'zpp_0000_0', 'z_0000_1', 'zp_0000_1', 'zpp_0000_1'],
          explain_rows = ['edge_0_0', 'edge_0_1', 'meridian_0_0', 'longitude_0_0'])

        The first row of the matrix means that the edge equation for
        edge 0 is

        .. math::

           {z_{0000,0}}^2 * z'_{0000,0} * z_{0000,1} * {z''_{0000,1}}^2 = 1.

        Similarly, the next row encodes the edge equation for the other edge
        and the next two rows encode peripheral equations.

        Following the SnapPy convention, a ``z`` denotes the cross ratio
        :math:`z` at the edge (0,1), a ``zp`` the cross ratio :math:`z'` at
        the edge (0,2) and a ``zpp`` the cross
        ratio :math:`z''` at the edge (1,2). The entire symbol ``z_xxxx_y`` then
        denotes the cross ratio belonging to the subsimplex at integral
        point ``xxxx`` (always ``0000`` for ``N = 2``) of the simplex ``y``.

        Note: the SnapPy convention is different from the paper
        mentioned above, e.g., compare
        kernel_code/edge_classes.c with Figure 3. We follow the SnapPy
        convention here so that all computations done in SnapPy are
        consistent.

        The explanations of the rows and columns can be obtained explicitly by:

        >>> M.gluing_equations_pgl(N = 3, equation_type = 'peripheral').explain_rows
        ['meridian_0_0', 'meridian_1_0', 'longitude_0_0', 'longitude_1_0']
        >>> M.gluing_equations_pgl(N = 2).explain_columns
        ['z_0000_0', 'zp_0000_0', 'zpp_0000_0', 'z_0000_1', 'zp_0000_1', 'zpp_0000_1']

        A subset of all gluing equations can be obtained by setting the
        ``equation_type``:

        * all gluing equations: ``all``
        * non-peripheral equations: ``non_peripheral``

          * edge gluing equations: ``edge``
          * face gluing equations: ``face``
          * internal gluing equations: ``internal``

        * cusp gluing equations: ``peripheral``

          * cusp gluing equations for meridians: ``meridian``
          * cusp gluing equations for longitudes: ``longitude``
        """

        cdef Integer_matrix_with_explanations c_matrix

        if N < 2 or N > 15:
            raise ValueError('N has to be 2...15')

        if equation_type not in ['all',
                                 'non_peripheral',
                                 'edge', 'face', 'internal',
                                 'peripheral',
                                 'longitude', 'meridian']:
            raise ValueError('Wrong equation_type')

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        equations = []
        explain_rows = []
        explain_cols = []

        if equation_type in [ 'all', 'non_peripheral', 'edge']:

            # Add edge equations
            get_edge_gluing_equations_pgl(self.c_triangulation,
                                          &c_matrix, N)
            eqns, r, explain_cols = convert_and_free_integer_matrix(c_matrix)
            equations += eqns
            explain_rows += r

        if equation_type in [ 'all', 'non_peripheral', 'face']:

            # Add face equations
            get_face_gluing_equations_pgl(self.c_triangulation,
                                          &c_matrix, N)
            eqns, r, explain_cols = convert_and_free_integer_matrix(c_matrix)
            equations += eqns
            explain_rows += r

        if equation_type in [ 'all', 'non_peripheral',  'internal']:

            # Add internal equations
            get_internal_gluing_equations_pgl(self.c_triangulation,
                                              &c_matrix, N)
            eqns, r, explain_cols = convert_and_free_integer_matrix(c_matrix)
            equations += eqns
            explain_rows += r

        if equation_type in ['all', 'peripheral', 'longitude', 'meridian']:

            # Add peripheral equations

            for i in range(self.num_cusps()):
                cusp_info = self.cusp_info(i)

                to_do = []
                # keep a todo list where we add (meridian, longitude)
                # pairs to process later

                if cusp_info.is_complete:

                    # Add meridians for complete cusps
                    if equation_type in [ 'meridian', 'peripheral', 'all']:
                        to_do += [ (1,0) ]
                        explain_rows += [
                            "meridian_%d_%d" % (j, i) for j in range(N-1) ]

                    # Add longitudes for complete cusps
                    if equation_type in [ 'longitude', 'peripheral', 'all']:
                        to_do += [ (0,1) ]
                        explain_rows += [
                            "longitude_%d_%d" % (j, i) for j in range(N-1) ]

                else:

                    # Add Dehn-filling for incomplete cusp
                    to_do += [ cusp_info.filling ]
                    explain_rows += [
                        "filling_%d_%d" % (j, i) for j in range(N-1) ]

                # process the todo list
                for (m, l) in to_do:

                    get_cusp_equations_pgl(self.c_triangulation,
                                           &c_matrix,
                                           N, i, m, l)

                    eqns, r, explain_cols = (
                        convert_and_free_integer_matrix(c_matrix))
                    equations += eqns

        if equations == []:
            # cover cases N = 2, 3 and equation_type = 'internal'
            return None

        return NeumannZagierTypeEquations(matrix(equations),
                                          explain_rows,
                                          explain_cols)

    def _ptolemy_equations_identified_face_classes(self):
        """
        This function returns an identification structure where s_f_t gets
        identified with -s_g_u if face f of tetrahedron t is glued to face g of
        tetrahedron u.

        We can represent a 2-cohomology class H^2(M,boundary M) by denoting by
        s_f_t the value the 2-cohomology class takes on the face f of
        tetrahedron t with the orientation being the one induced from the
        orientation of the tetrahedron.
        Because a face class of the triangulation has two representatives
        (tet_index, face_index) and the gluing is orientation-reversing on the
        face, one s will be the negative of another s.
        """

        cdef Identification_of_variables c_vars

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        number_the_tetrahedra(self.c_triangulation)

        get_ptolemy_equations_identified_face_classes(
            self.c_triangulation, &c_vars)

        return convert_and_free_identification_of_variables(c_vars)

    def _ptolemy_equations_identified_coordinates(self, N,
                                                  obstruction_class = None):

        """
        Ptolemy coordinates that need to be identified for the given
        triangulation when computing pSL(N,C) representations.
        """

        cdef Identification_of_variables c_vars
        cdef int *c_obstruction_class = NULL
        cdef int T

        if N < 2 or N > 15:
            raise ValueError('N has to be 2...15')

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        if obstruction_class:
            T = get_num_tetrahedra(self.c_triangulation)
            if 2 * T != len(obstruction_class):
                raise ValueError('Obstruction class has wrong length')
            c_obstruction_class = <int *>malloc(2*T*sizeof(int))
            for i, c in enumerate(obstruction_class):
                c_obstruction_class[i] = c

        number_the_tetrahedra(self.c_triangulation)

        get_ptolemy_equations_identified_coordinates(
            self.c_triangulation, &c_vars, N, c_obstruction_class)

        free(c_obstruction_class)

        return convert_and_free_identification_of_variables(c_vars)

    def _ptolemy_equations_action_by_decoration_change(self, int N):
        """
        We can change a decoration by multiplying a coset of a cusp by a
        diagonal matrix. Let a diagonal matrix SL(n,C) with diagonal
        entries 1 1 ... z 1 ... 1 1/z (z at position j) act on cusp i. It
        changes some Ptolemy coordinate c_p_t by some power z^n.
        This is expressed in the following matrix as the entry in row
        labelled c_p_t and the column labelled diagonal_entry_j_on_cusp_i.
        """

        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        number_the_tetrahedra(self.c_triangulation)

        get_ptolemy_equations_action_by_decoration_change(
            self.c_triangulation,N, &c_matrix)
        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def _ptolemy_equations_boundary_map_3(self):
        """
        Boundary map C_3 -> C_2 in cellular homology represented as matrix

        The following map represents the boundary map in the cellular chain
        complex when representing a linear map as a matrix m acting on a column
        vector v by left-multiplication m * v. With right-multiplication acting
        on row vectors, the matrix represents the coboundary map in the cochain
        complex.

        The basis for C_3 are just the oriented tetrahedra of the triangulation.
        The basis for C_2 are the face classes, see
        _ptolemy_equations_identified_face_classes.
        """

        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        number_the_tetrahedra(self.c_triangulation)

        get_ptolemy_equations_boundary_map_3(self.c_triangulation, &c_matrix)

        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def _ptolemy_equations_boundary_map_2(self):
        """
        Boundary map C_2 -> C_1 in cellular homology represented as matrix.

        Also see _ptolemy_equations_boundary_map_3.
        """

        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        number_the_tetrahedra(self.c_triangulation)

        get_ptolemy_equations_boundary_map_2(self.c_triangulation, &c_matrix)

        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def _ptolemy_equations_boundary_map_1(self):
        """
        Boundary map C_1 -> C_0 in cellular homology represented as matrix.
        This will compute the homology of the cell complex obtained when
        gluing together the tetrahedra and not of the cusped manifold.

        Also see _ptolemy_equations_boundary_map_3.
        """

        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        number_the_tetrahedra(self.c_triangulation)

        get_ptolemy_equations_boundary_map_1(self.c_triangulation, &c_matrix)

        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def ptolemy_obstruction_classes(self):

        """
        Returns the obstruction classes needed to compute
        pSL(N,C) = SL(N,C)/{+1,-1} representations for even N, i.e., it
        returns a list with a representative cocycle for each class in
        H^2(M, boundary M; Z/2). The first element in the list is always
        representing the trivial obstruction class.

        For example, 4_1 has two obstruction classes:

        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_obstruction_classes()
        >>> len(c)
        2

        The primary use of these obstruction classes is to construct
        the Ptolemy variety as described in Definition 1.7 of
        Stavros Garoufalidis, Dylan Thurston, Christian K. Zickert:
        "The Complex Volume of SL(n,C)-Representations of 3-Manifolds"
        (http://arxiv.org/abs/1111.2828).

        For example, to construct the Ptolemy variety for
        PSL(2,C)-representations of 4_1 that do not lift to boundary-parabolic
        SL(2,C)-representations, use:

        >>> p = M.ptolemy_variety(N = 2, obstruction_class = c[1])

        Or the following short-cut:

        >>> p = M.ptolemy_variety(2, obstruction_class = 1)

        Note that this obstruction class only makes sense for even N:

        >>> p = M.ptolemy_variety(3, obstruction_class = c[1])
        Traceback (most recent call last):
        ...
        AssertionError: PtolemyObstructionClass only makes sense for even N, try PtolemyGeneralizedObstructionClass

        To obtain PGL(N,C)-representations for N > 2, use the generalized
        obstruction class:

        >>> c = M.ptolemy_generalized_obstruction_classes(3)
        >>> p = M.ptolemy_variety(3, obstruction_class = c[1])

        The original obstruction class encodes a representing cocycle in Z/2 as follows:

        >>> c = M.ptolemy_obstruction_classes()
        >>> c[1]
        PtolemyObstructionClass(s_0_0 + 1, s_1_0 - 1, s_2_0 - 1, s_3_0 + 1, s_0_0 - s_0_1, s_1_0 - s_3_1, s_2_0 - s_2_1, s_3_0 - s_1_1)

        This means that the cocycle to represent this obstruction class in Z/2
        takes value 1 in Z/2 on face 0 of tetrahedra 0 (because s_0_0 = -1)
        and value 0 in Z/2 on face 1 of tetrahedra 0 (because s_1_0 = +1).

        Face 3 of tetrahedra 0 and face 1 of tetrahedra 1 are identified,
        hence the cocycle takes the same value on those two faces (s_3_0 = s_1_1).

        """

        return ptolemyManifoldMethods.get_ptolemy_obstruction_classes(self)

    def ptolemy_generalized_obstruction_classes(self, N):

        """
        Returns the obstruction classes needed to compute
        PGL(N,C)-representations for any N, i.e., it returns a list with
        a representative cocycle for each element in
        H^2(M, boundary M; Z/N) / (Z/N)^* where (Z/N)^* are the units in Z/N.
        The first element in the list always corresponds to the trivial
        obstruction class.
        The generalized ptolemy obstruction classes are thus a generalization
        of the ptolemy obstruction classes that allow to find all
        boundary-unipotent
        PGL(N,C)-representations including those that do not lift to
        boundary-unipotent SL(N,C)-representations for N odd or
        SL(N,C)/{+1,-1}-representations for N even.

        For example, 4_1 has three obstruction classes up to equivalence:

        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_generalized_obstruction_classes(4)
        >>> len(c)
        3

        For 4_1, we only get three obstruction classes even though we have
        H^2(M, boundary M; Z/4) = Z/4 because the two obstruction classes
        1 in Z/4 and -1 in Z/4 are related by a unit and thus give
        isomorphic Ptolemy varieties.

        The primary use of an obstruction class sigma is to construct the
        Ptolemy variety of sigma. This variety computes boundary-unipotent
        PGL(N,C)-representations whose obstruction class to a
        boundary-unipotent lift to SL(N,C) is sigma.

        For example for 4_1, there are 2 obstruction classes for N = 3:

        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_generalized_obstruction_classes(3)
        >>> len(c)
        2

        The Ptolemy variety parametrizing boundary-unipotent
        SL(3,C)-representations of 4_1 is obtained by

        >>> p = M.ptolemy_variety(N = 3, obstruction_class = c[0])

        and the Ptolemy variety parametrizing boundary-unipotent
        PSL(3,C)-representations of 4_1 that do not lift to
        boundary-unipotent SL(3,C)-representations is obtained by

        >>> p = M.ptolemy_variety(N = 3, obstruction_class = c[1])

        The cocycle representing the non-trivial obstruction class looks as
        follows:

        >>> c[1]
        PtolemyGeneralizedObstructionClass([2, 0, 0, 1])

        This means that the cocycle takes the value -1 in Z/3 on the first face
        class and 1 on the fourth face class but zero on every other of the
        four face classes.
        """
        return (
            ptolemyManifoldMethods.get_generalized_ptolemy_obstruction_classes(
                self, N))

    def ptolemy_variety(self, N, obstruction_class = None,
                        simplify = True, eliminate_fixed_ptolemys = False):

        """
        M.ptolemy_variety(N, obstruction_class = None, simplify = True, eliminate_fixed_ptolemys = False)

        Returns a Ptolemy variety as described in

        * Stavros Garoufalidis, Dyland Thurston, Christian K. Zickert:
          "The Complex Volume of SL(n,C)-Representations of 3-Manifolds"
          (http://arxiv.org/abs/1111.2828)
        * Stavros Garoufalidis, Matthias Goerner, Christian K. Zickert:
          "Gluing Equations for PGL(n,C)-Representations of 3-Manifolds "
          (http://arxiv.org/abs/1207.6711)

        The variety can be exported to magma or sage and solved there. The
        solutions can be processed to compute invariants. The method can also
        be used to automatically look up precomputed solutions from the
        database at http://ptolemy.unhyperbolic.org/data .

        Example for m011 and PSL(2,C)-representations:

        >>> M = Manifold("m011")

        Obtain all Ptolemy varieties for PSL(2,C)-representations:

        >>> p = M.ptolemy_variety(2, obstruction_class = 'all')

        There are two Ptolemy varieties for the two obstruction classes:

        >>> len(p)
        2

        Retrieve the solutions from the database

        >>> sols = p.retrieve_solutions() #doctest: +SKIP

        Compute the solutions using magma (default in SnapPy)

        >>> sols = p.compute_solutions(engine = 'magma') #doctest: +SKIP

        Compute the solutions using singular (default in sage)

        >>> sols = p.compute_solutions(engine = 'sage') #doctest: +SKIP

        Note that magma is significantly faster.

        Compute all resulting complex volumes

        >>> cvols = sols.complex_volume_numerical() #doctest: +SKIP
        >>> cvols  #doctest: +SKIP
        [[[-4.29405713186238 E-16 + 0.725471193740844*I,
           -0.942707362776931 + 0.459731436553693*I,
           0.942707362776931 + 0.459731436553693*I]],
         [[3.94159248086745 E-15 + 0.312682687518267*I,
           4.64549527022581 E-15 + 0.680993020093457*I,
           -2.78183391239608 - 0.496837853805869*I,
           2.78183391239608 - 0.496837853805869*I]]]

        Show complex volumes as a non-nested list:

        >>> cvols.flatten(depth=2) #doctest: +SKIP
        [-4.29405713186238 E-16 + 0.725471193740844*I,
         -0.942707362776931 + 0.459731436553693*I,
         0.942707362776931 + 0.459731436553693*I,
         3.94159248086745 E-15 + 0.312682687518267*I,
         4.64549527022581 E-15 + 0.680993020093457*I,
         -2.78183391239608 - 0.496837853805869*I,
         2.78183391239608 - 0.496837853805869*I]

        For more examples, go to http://ptolemy.unhyperbolic.org/

        === Optional Arguments ===

        obstruction_class --- class from Definition 1.7 of (1).
        None for trivial class or a value returned from ptolemy_obstruction_classes.
        Short cuts: obstruction_class = 'all' returns a list of Ptolemy varieties
        for each obstruction. For easier iteration, can set obstruction_class to
        an integer.

        simplify --- boolean to indicate whether to simplify the equations which
        significantly reduces the number of variables.
        Simplifying means that several identified Ptolemy coordinates x = y = z = ...
        are eliminated instead of adding relations x - y = 0, y - z = 0, ...

        eliminate_fixed_ptolemys --- boolean to indicate whether to eliminate
        the Ptolemy coordinates that are set to 1 for fixing the decoration.
        Even though this simplifies the resulting representation, setting it to
        True can cause magma to run longer when finding a Groebner basis.

        === Examples for 4_1 ===

        >>> M = Manifold("4_1")

        Get the varieties for all obstruction classes at once (use
        help(varieties[0]) for more information):

        >>> varieties = M.ptolemy_variety(2, obstruction_class = "all")

        Print the variety as an ideal (sage object) for the non-trivial class:

        >>> varieties[1].ideal    #doctest: +SKIP
        Ideal (-c_0011_0^2 + c_0011_0*c_0101_0 + c_0101_0^2, -c_0011_0^2 - c_0011_0*c_0101_0 + c_0101_0^2, c_0011_0 - 1) of Multivariate Polynomial Ring in c_0011_0, c_0101_0 over Rational Field

        Print the equations of the variety for the non-trivial class:

        >>> for eqn in varieties[1].equations:
        ...     print(eqn)          #doctest: +NORMALIZE_WHITESPACE
             - c_0011_0 * c_0101_0 + c_0011_0^2 + c_0101_0^2
             c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
             - 1 + c_0011_0

        Generate a magma file to compute Primary Decomposition for N = 3:

        >>> p = M.ptolemy_variety(3)
        >>> s = p.to_magma()
        >>> print(s.split("ring and ideal")[1].strip())     #doctest: +ELLIPSIS
        R<c_0012_0, c_0012_1, c_0102_0, c_0111_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0> := PolynomialRing(RationalField(), 8, "grevlex");
        MyIdeal := ideal<R |
                  c_0012_0 * c_1101_0 + c_0102_0 * c_0111_0 - c_0102_0 * c_1011_0,
            ...

        === If you have a magma installation ===

        Call p.compute_solutions() to automatically call magma on the above output
        and produce exact solutions!!!

        >>> try:
        ...     sols = p.compute_solutions()
        ... except:
        ...     # magma failed, use precomputed_solutions
        ...     sols = None

        Check solutions against manifold
        >>> if sols:
        ...     dummy = sols.check_against_manifold()

        === If you do not have a magma installation ===

        Load a precomputed example from magma which is provided with the package:

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
        >>> print(_magma_output_for_4_1__sl3)      #doctest: +ELLIPSIS
        <BLANKLINE>
        ==TRIANGULATION=BEGINS==
        % Triangulation
        4_1
        ...

        Parse the file and produce solutions:

        >>> sols = solutions_from_magma(_magma_output_for_4_1__sl3)

        >>> dummy = sols.check_against_manifold()

        === Continue here whether you have or do not have magma ===

        Pick the first solution of the three different solutions (up to Galois
        conjugates):

        >>> len(sols)
        3
        >>> solution = sols[0]

        Read the exact value for c_1020_0 (help(solution) for more information
        on how to compute cross ratios, volumes and other invariants):

        >>> solution['c_1020_0']
        Mod(-1/2*x - 3/2, x^2 + 3*x + 4)

        Example of simplified vs non-simplified variety for N = 4:

        >>> simplified = M.ptolemy_variety(4, obstruction_class = 1)
        >>> full = M.ptolemy_variety(4, obstruction_class = 1, simplify = False)
        >>> len(simplified.variables), len(full.variables)
        (21, 63)
        >>> len(simplified.equations), len(full.equations)
        (24, 72)
        """

        return ptolemyManifoldMethods.get_ptolemy_variety(
            self, N, obstruction_class,
            simplify = simplify,
            eliminate_fixed_ptolemys = eliminate_fixed_ptolemys)

    def gluing_equations(self,form='log'):
        """
        In the default mode, this function returns a matrix with rows
        of the form

                  a b c  d e f  ...

        which means

            a*log(z0) + b*log(1/(1-z0)) + c*log((z0-1)/z0) + d*log(z1) +... = 2 pi i

        for an edge equation, and (same) = 0 for a cusp equation.
        Here, the cusp equations come at the bottom of the matrix, and
        are listed in the form: meridian of cusp 0, longitude of cusp
        0, meridian of cusp 1, longitude of cusp 1,...

        In terms of the tetrahedra, a is the invariant of the edge
        (2,3), b the invariant of the edge (0,2) and c is the
        invariant of the edge (1,2).  See kernel_code/edge_classes.c
        for a detailed account of the convention used.

        If the optional argument form='rect' is given, then this
        function returns a list of tuples of the form:

           ( [a0, a1,..,a_n], [b_0, b_1,...,b_n], c)

        where this corresponds to the equation

           z0^a0 (1 - z0)^b0 z1^a1(1 - z1)^b1 ...  = c

        where c = 1 or -1.

        >>> M = Triangulation('m004(2,3)')
        >>> M.gluing_equations()
        [ 2  1  0  1  0  2]
        [ 0  1  2  1  2  0]
        [ 2  0  0  0 -8  6]
        >>> M.gluing_equations(form='rect')
        [([2, -1], [-1, 2], 1), ([-2, 1], [1, -2], 1), ([2, -6], [0, 14], 1)]
        """

        cdef int **c_eqns
        cdef int num_rows, num_cols
        cdef int* eqn

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        c_eqns = get_gluing_equations(self.c_triangulation,
                                      &num_rows, &num_cols)
        eqns = [ [ c_eqns[i][j] for j in range(num_cols) ]
                 for i in range(num_rows) ]
        free_gluing_equations(c_eqns, num_rows)

        for i in range(self.num_cusps()):
            cusp_info = self.cusp_info(i)
            if cusp_info.is_complete:
                to_do = [(1,0), (0,1)]
            else:
                to_do = [cusp_info.filling]
            for (m, l) in to_do:
                eqn = get_cusp_equation(self.c_triangulation,
                                        i, int(m), int(l), &num_rows)
                eqns.append([eqn[j] for j in range(num_rows)])
                free_cusp_equation(eqn)

        if form == 'log':
            return matrix(eqns)

        if form != 'rect':
            raise ValueError("Equations are available in 'log' and "
                             "'rect' forms only.")
        rect = []
        for row in eqns:
            n = self.num_tetrahedra()
            a, b = [0,]*n, [0,]*n
            c = 1
            for j in range(n):
                r = row[3*j + 2]
                a[j] = row[3*j] - r
                b[j] = -row[3*j + 1] + r
                c *= -1 if r % 2 else 1
            rect.append( (a, b, c) )
        return rect

    cdef big_homology(self):
        """
        Returns an AbelianGroup representing the first integral
        homology group of the underlying (Dehn filled) manifold.
        Preliminary simplification is done with arbitrary precision
        integers.  Smith form is then computed with PARI.

        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z
        """

        relation_matrix = self.c_presentation_matrix()
        return AbelianGroup(relation_matrix.simplified_matrix())

    cdef c_presentation_matrix(self):
        """
        See _presentation_matrix.
        """

        if not all_Dehn_coefficients_are_integers(self.c_triangulation):
            raise ValueError('All Dehn filling coefficients must be integers')
        choose_generators(self.c_triangulation, 0, 0)
        num_generators = self.c_triangulation.num_generators
        cdef EdgeClass* edge
        cdef PositionedTet ptet, ptet0
        cdef GeneratorStatus status
        cdef c_Tetrahedron* tet
        cdef VertexIndex vertex
        cdef FaceIndex side
        cdef Orientation orientation
        cdef int row, column, num_edges, m, l
        relation_matrix = PresentationMatrix(0, num_generators)
        edge = self.c_triangulation.edge_list_begin.next
        # find the edge relations
        while edge != &(self.c_triangulation.edge_list_end):
            row = relation_matrix.rows
            relation_matrix.add_rows(1)
            set_left_edge(edge, &ptet0)
            ptet = ptet0
            while True:
                column = ptet.tet.generator_index[ptet.near_face]
                status = ptet.tet.generator_status[ptet.near_face]
                if status == outbound_generator:
                    relation_matrix[row, column] += 1
                elif status == inbound_generator:
                    relation_matrix[row, column] -= 1
                elif status != not_a_generator:
                    raise RuntimeError('Invalid generator status')
                veer_left(&ptet)
                if same_positioned_tet(&ptet, &ptet0):
                    break
            row += 1
            edge = edge.next
        # find the cusp relations
        num_edges = relation_matrix.rows
        relation_matrix.add_rows(self.num_cusps())
        tet = self.c_triangulation.tet_list_begin.next
        while tet != &(self.c_triangulation.tet_list_end):
            for vertex in range(4):
                if tet.cusp[vertex].is_complete:
                    continue
                for side in range(4):
                    if side == vertex:
                        continue
                    if tet.generator_status[side] != inbound_generator:
                        continue
                    for orientation in (right_handed, left_handed):
                        row = num_edges + tet.cusp[vertex].index
                        column = tet.generator_index[side]
                        m = <int>tet.cusp[vertex].m
                        l = <int>tet.cusp[vertex].l
                        relation_matrix[row, column] += (
                            m * tet.curve[0][orientation][vertex][side]
                            + l * tet.curve[1][orientation][vertex][side])
            tet = tet.next

        return relation_matrix

    def _presentation_matrix(self):
        """
        Computes the presentation matrix for the first integral
        homology group of the underlying (Dehn filled) manifold.
        """
        return self.c_presentation_matrix()

    cdef csmall_homology(self):
        """
        Returns an AbelianGroup representing the first integral
        homology group of the underlying (Dehn filled) manifold.
        Preliminary simplification is done with 64 bit integers.
        Smith form is then computed with PARI.

        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z
        """
        cdef RelationMatrix R
        cdef int n

        if self.c_triangulation is NULL:
            return AbelianGroup()
        homology_presentation(self.c_triangulation, &R)
        relations = []
        if R.relations != NULL:
            if R.num_rows == 0:
                relations = [0,] * R.num_columns
            else:
                for m from 0 <= m < R.num_rows:
                    row = []
                    for n from 0 <= n < R.num_columns:
                        row.append(R.relations[m][n])
                    relations.append(row)
                free_relations(&R)
        else:
            raise RuntimeError("The SnapPea kernel couldn't compute "
                               "the homology presentation matrix")
        return AbelianGroup(relations)

    def homology(self) -> AbelianGroup:
        """
        Returns an :class:`AbelianGroup` representing the first integral
        homology group of the underlying (Dehn filled) manifold.

        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z

        """
        try:
            return self._cache.lookup('homology')
        except KeyError:
            pass

        cdef c_AbelianGroup *H
        cdef int n

        if self.c_triangulation is NULL:
            return AbelianGroup()
        H = homology(self.c_triangulation)
        if H != NULL:
            coefficient_list = []
            compress_abelian_group(H)
            for n from 0 <= n < H.num_torsion_coefficients:
                coefficient_list.append(H.torsion_coefficients[n])
            free_abelian_group(H)
            result = AbelianGroup(elementary_divisors=coefficient_list)
        else:
            try:
                result = self.csmall_homology()
            except RuntimeError:
                result = self.big_homology()
        return self._cache.save(result, 'homology')

    def fundamental_group(
            self,
            simplify_presentation : bool = True,
            fillings_may_affect_generators : bool = True,
            minimize_number_of_generators : bool = True,
            try_hard_to_shorten_relators : bool = True
            ) -> FundamentalGroup:
        """
        Returns a :class:`FundamentalGroup` object representing the fundamental
        group of the manifold.  If integer Dehn surgery parameters
        have been set, then the corresponding peripheral elements are
        killed.

        >>> M = Triangulation('m004')
        >>> G = M.fundamental_group()
        >>> G
        Generators:
           a,b
        Relators:
           aaabABBAb
        >>> G.peripheral_curves()
        [('ab', 'aBAbABab')]

        There are four optional arguments all of which default to True:

        - simplify_presentation
        - fillings_may_affect_generators
        - minimize_number_of_generators
        - try_hard_to_shorten_relators

        >>> M.fundamental_group(False, False, False)
        Generators:
           a,b,c
        Relators:
           CbAcB
           BacA
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        args = (simplify_presentation, fillings_may_affect_generators,
                minimize_number_of_generators, try_hard_to_shorten_relators)
        try:
            return self._cache.lookup('fundamental_group', *args)
        except KeyError:
            pass

        return self._cache.save(FundamentalGroup(self, *args),
                                'fundamental_group', *args)

    def cover(self, permutation_rep) -> 'snappy.Triangulation':
        """
        Returns a :class:`Triangulation` representing the finite cover specified
        by a transitive permutation representation.  The representation is
        specified by a list of permutations, one for each generator of the
        simplified presentation of the fundamental group.  Each permutation is
        specified as a list ``P`` such that ``set(P) == set(range(d))`` where
        ``d`` is the degree of the cover.

        >>> M = Triangulation('m004')
        >>> N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])
        >>> N0.homology()
        Z + Z + Z
        >>> N0.cover_info()['type']
        'irregular'
        >>> N0.cover_info()['base']
        'm004'
        >>> N0.cover_info()['degree']
        5

        Within Sage the permutations can also be of type
        ``PermutationGroupElement``, in which case they act on the set
        ``range(1, d + 1)``.  Or, you can specify a GAP or Magma subgroup
        of the fundamental group.  For more examples, see the docstring for
        :meth:`Manifold.cover`.
        """
        cdef RepresentationIntoSn* c_representation
        cdef c_Triangulation* c_triangulation
        cdef Triangulation cover

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        # For Sage, we need to check if we have been given some
        # alternate inputs

        if _within_sage:
            if is_GapElement(permutation_rep):
                if permutation_rep.IsSubgroupFpGroup():
                    GG = gap(self.fundamental_group())
                    coset_action = GG.FactorCosetAction(permutation_rep)
                    gap_image_gens = coset_action.Image().GeneratorsOfGroup()
                    Q = PermutationGroup(gap_image_gens)
                    return self.cover([Q(g) for g in gap_image_gens ])
                elif (permutation_rep.IsToPermGroupHomomorphismByImages() or
                      permutation_rep.IsCompositionMappingRep()):
                    f = permutation_rep
                    return self.cover(f.PreImage(f.Image().Stabilizer(1)))
                else:
                    raise ValueError('Cover specified by unrecognized GAP type')

            elif is_MagmaElement(permutation_rep):
                input_type = repr(permutation_rep.Type())
                if input_type == 'GrpFP':
                    GG = magma(self.fundamental_group())
                    f = GG.CosetAction(permutation_rep)
                elif input_type == 'HomGrp':
                    f = permutation_rep
                    if not repr(f.Image().Type()) == 'GrpPerm':
                        raise TypeError('The homomorphism image is not '
                                        'a permutation group.')
                else:
                    raise TypeError('That Magma type not recognized.')

                magma.eval("""\
                     FormatHomForSnapPea := function(f)
                         subone := function(L)   return [x - 1 : x in L]; end function;
                         return [subone(Eltseq(f(g))) : g in Generators(Domain(f))];
                       end function;""")
                permutation_rep = f.FormatHomForSnapPea().sage()

            # Not a useful GAP or MAGMA object, so let's try.
            elif all(isinstance(p, PermutationGroupElement) for p in permutation_rep):
                permutation_rep = [[x - 1 for x in perm.domain()]
                                   for perm in permutation_rep]

        G = self.fundamental_group()
        c_representation = self.build_rep_into_Sn(permutation_rep)
        degree = len(permutation_rep[0]) if len(permutation_rep) > 0 else 1

        # The next call has the effect of initializing aspects of
        # self.c_triangulation that are needed to build the cover.  It is *not*
        # redundant with the preceeding call to self.fundamental_group()
        # because the resulting FundamentalGroup object calls the kernel
        # function "fundamental_group" on a copy of self.c_triangulation.
        free_group_presentation(compute_unsimplified_presentation(self.c_triangulation))

        c_triangulation = construct_cover(self.c_triangulation,
                                          c_representation,
                                          <int>degree)
        cover = self.__class__('empty')
        cover.set_c_triangulation(c_triangulation)
        cover._cover_info = info = {
            'base'   : self.name(),
            'type'   : cover_types[c_representation.covering_type],
            'degree' : degree
            }
        cover.set_name(self.name() +'~' + info['type'][:3])
        free_representation(c_representation,
                            G.num_original_generators(),
                            self.num_cusps())
        return cover

    def covers(
            self,
            degree,
            method : typing.Optional[str] = None,
            cover_type : str ='all'
        ) -> 'list[snappy.Triangulation]':
        """
        Returns a list of Triangulations corresponding to all of the
        finite covers of the given degree.  The default method is
        'low_index' for general covers and 'snappea' for cyclic
        covers.  The former uses Sim's algorithm while the latter
        uses the original Snappea algorithm.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        >>> M = Triangulation('m003')
        >>> covers = M.covers(4)
        >>> sorted(N.homology() for N in covers)
        [Z/3 + Z/15 + Z, Z/5 + Z + Z]

        It is faster to look just at cyclic covers.

        >>> covers = M.covers(4, cover_type='cyclic')
        >>> [(N, N.homology()) for N in covers]
        [(m003~cyc~0(0,0), Z/3 + Z/15 + Z)]

        Here we check that we get the same number of covers with the
        'snappea' and 'low_index' methods.

        >>> M = Triangulation('m125')
        >>> len(M.covers(5))
        19
        >>> len(M.covers(5, method='snappea'))
        19

        If you are using Sage, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional
        argument method = 'gap' If you have Magma installed, you can
        used it to do the heavy lifting by specifying method='magma'.
        """
        if degree < 1:
            raise ValueError('Cover degree should be at least 1')

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if cover_type not in ('cyclic', 'all'):
            raise ValueError("Supported cover_types are 'all' "
                             "and 'cyclic'.")
        if method is None:
            method = 'low_index'
        else:
            method = method.lower()

        if cover_type == 'cyclic' or method == 'snappea':
            return self._covers_snappea(degree, cover_type)

        if method == 'low_index':
            return self._covers_low_index(degree)
        if method == 'gap':
            return self._covers_gap(degree)
        if method == 'magma':
            return self._covers_magma(degree)

        raise ValueError("Supported methods are 'low_index', 'gap', 'magma' "
                         "and 'snappea'")

    def _covers_low_index(self, degree, num_threads=0):
        """
        Compute all covers using low_index.
        """
        G = self.fundamental_group()

        relators = sorted(G.relators(as_int_list=True), key=len)

        short_relators, long_relators = [], relators
        strategy = ''

        if relators:
            if len(relators[0]) <= 2*degree:
                short_relators = [rel for rel in relators if len(rel) <= 2*degree]
                long_relators = relators[len(short_relators):]
                strategy = 'spin_short'
            elif len(relators[0]) <= 3*degree:
                short_relators = [rel for rel in relators if len(rel) <= 3*degree]
                long_relators = relators[len(short_relators):]

        reps = low_index.permutation_reps(G.num_generators(),
                                          short_relators,
                                          long_relators,
                                          degree,
                                          strategy=strategy,
                                          num_threads=num_threads)

        def index(subgroup):
            return 1 if len(subgroup) == 0 else len(subgroup[0])

        return [self.cover(H) for H in reps if index(H) == degree]

    def _covers_gap(self, degree):
        """
        Compute all covers using gap.
        """
        if not _within_sage:
            raise SageNotAvailable('the "gap" method for covers requires Sage')

        G = gap(self.fundamental_group())
        return [self.cover(H)
                for H in G.LowIndexSubgroupsFpGroup(degree)
                if G.Index(H) == degree]

    def _covers_magma(self, degree):
        """
        Compute all covers using magma.
        """
        if not _within_sage:
            raise SageNotAvailable('the "magma" method for covers requires Sage')

        G = magma(self.fundamental_group())
        return [self.cover(H)
                for H in G.LowIndexSubgroups('<%d, %d>' %
                                             (degree, degree))]

    def _covers_snappea(self, degree, cover_type):
        """
        Compute cyclic or all covers using SnapPea kernel.
        """
        cdef RepresentationList* reps
        cdef RepresentationIntoSn* rep
        cdef c_Triangulation* cover
        cdef Triangulation T
        cdef PermutationSubgroup c_cover_type

        if cover_type == 'cyclic':
            c_cover_type = permutation_subgroup_Zn
        else:
            c_cover_type = permutation_subgroup_Sn

        reps = find_representations(
            self.c_triangulation, degree, c_cover_type)

        covers = []
        rep = reps.list
        cover_count = 0
        while rep != NULL:
            cover = construct_cover(self.c_triangulation,
                                    rep,
                                    reps.num_sheets)
            T = self.__class__('empty')
            T.set_c_triangulation(cover)
            T._cover_info = info = {
                'base'   : self.name(),
                'type'   : cover_types[rep.covering_type],
                'degree' : degree
                }
            T.set_name(info['base'] + "~" + info['type'][:3] + '~%d' %
                       cover_count)
            covers.append(T)
            cover_count += 1
            rep = rep.next

        free_representation_list(reps)
        return covers

    cdef RepresentationIntoSn *build_rep_into_Sn(self, perm_list) except ? NULL:
        """
        Build a SnapPea RepresentationIntoSn from a list of
        permutations, one for each generator of the simplified
        fundamental group.  A permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.  The representation constructed here is given in terms
        of the geometric generators, for use in constructing a covering
        space.
        """
        cdef c_Triangulation* c_triangulation
        cdef RepresentationIntoSn* c_representation
        cdef RepresentationIntoSn* c_repn_in_original_gens = NULL
        cdef int i, j
        cdef int** c_original_generators
        cdef int** c_relators
        cdef int** c_meridians
        cdef int** c_longitudes

        if len(perm_list) == 0:  # implies trivial presentation
            degree = 1
        else:
            degree = len(perm_list[0])

        # Sanity check
        S = set(range(degree))
        for permutation in perm_list:
            if set(permutation) != S:
                raise ValueError('The permutation list is invalid.')

        # Initialize
        num_cusps = self.num_cusps()
        c_triangulation = self.c_triangulation
        G = self.fundamental_group()
        num_generators = G.num_generators()
        if len(perm_list) != num_generators:
            raise ValueError('Number of permutations is not the same '
                             'as the number of generators')

        relators = G.relators(as_int_list=True)
        num_relators = len(relators)
        orig_gens = G.original_generators(as_int_list=True)
        num_orig_gens = len(orig_gens)
        peripheral_curves = G.peripheral_curves(as_int_list=True)

        # Allocate a whole bunch of memory, SnapPea and malloc.
        c_representation = initialize_new_representation(
            num_orig_gens,
            <int>degree,
            num_cusps)
        for i from 0 <= i < num_generators:
            for j from 0 <= j < degree:
                c_representation.image[i][j] = perm_list[i][j]
        c_original_generators = <int**>malloc(num_orig_gens*sizeof(int*))
        for i from  0 <= i < num_orig_gens:
            c_original_generators[i] = c_word_from_list(orig_gens[i])
        c_relators = <int**>malloc(num_relators*sizeof(int*))
        for i from  0 <= i < num_relators:
            c_relators[i] = c_word_from_list(relators[i])
        c_meridians = <int**>malloc(num_cusps*sizeof(int*))
        c_longitudes = <int**>malloc(num_cusps*sizeof(int*))
        for i from 0 <= i < num_cusps:
            meridian, longitude = peripheral_curves[i]
            c_meridians[i] = c_word_from_list(meridian)
            c_longitudes[i] = c_word_from_list(longitude)
        # Whew!

        failed = False
        if (candidateSn_is_valid(c_representation.image,
                                 <int>degree, c_relators, num_relators) and
            candidateSn_is_transitive(c_representation.image,
                                      num_generators, <int>degree) ):
            c_repn_in_original_gens = convert_candidateSn_to_original_generators(
                c_representation.image,
                <int>degree,
                num_orig_gens,
                c_original_generators,
                c_triangulation,
                c_meridians,
                c_longitudes)
        else:
            message = 'Invalid permutation data.'
            failed = True
        if c_repn_in_original_gens == NULL:
            message = 'Failed to construct permutation representation.'
            failed = True

        # Now free all that memory
        for i from 0 <= i < num_cusps:
            free(c_meridians[i])
            free(c_longitudes[i])
        free(c_meridians)
        free(c_longitudes)
        for i from 0 <= i < num_relators:
            free(c_relators[i])
        free(c_relators)
        for i from 0 <= i < num_orig_gens:
            free(c_original_generators[i])
        free(c_original_generators)
        free_representation(c_representation, num_generators, num_cusps)
        # Free at last!

        if failed:
            raise RuntimeError(message)
        return c_repn_in_original_gens

    def set_peripheral_curves(self, peripheral_data,
                              which_cusp=None, return_matrices=False):
        """
        Each cusp has a preferred marking. In the case of a torus
        cusp, this is pair of essential simple curves meeting in one
        point; equivalently, a basis of the first homology of the
        boundary torus. These curves are called the meridian and the
        longitude.

        This method changes these markings in various ways.  In many
        cases, if the flag return_matrices is True then it returns
        a list of change-of-basis matrices is returned, one per
        cusp, which will restore the original markings if passed
        as peripheral_data.
        """
        cdef int a,b,c,d
        cdef MatrixInt22 *matrices
        cdef c_FuncResult result

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty')

        if which_cusp is not None:
            which_cusp = valid_index(
                which_cusp, self.num_cusps(),
                'The specified cusp (%s) does not exist.')

        self._cache.clear(message='set_peripheral_curves')
        if peripheral_data == 'fillings':
            if which_cusp is not None:
                raise ValueError("You must apply 'fillings' to all "
                                 "of the cusps.")
            install_current_curve_bases(self.c_triangulation)
            return
        elif peripheral_data == 'combinatorial':
            if return_matrices:
                matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                                 sizeof(MatrixInt22))
                install_combinatorial_bases(self.c_triangulation, matrices)
                cobs = []
                for n in range(self.num_cusps()):
                    cobs.append([ [ matrices[n][0][0], matrices[n][0][1] ],
                                  [ matrices[n][1][0], matrices[n][1][1] ] ])
                free(matrices)
                return cobs
            else:
                peripheral_curves(self.c_triangulation)
        elif which_cusp is not None:
            meridian, longitude = peripheral_data
            a, b = meridian
            c, d = longitude
            if a*d - b*c != 1 and (self.is_orientable() or a*d - b*c != -1):
                if self.is_orientable():
                    raise ValueError('The data provided does not give a '
                                     '(positively oriented) basis.')
                else:
                    raise ValueError('The data provided does not give a basis.')

            if 'Klein' in self.cusp_info(which_cusp)['topology']:
                if b != 0 or c != 0:
                    raise ValueError('The data provided are invalid for Klein '
                                     'bottle cusps.')

            matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                             sizeof(MatrixInt22))

            for n in range(self.num_cusps()):
                for i,j in [(0,0),(0,1),(1,0),(1,1)]:
                    matrices[n][i][j] = 1 if i == j else 0

            matrices[which_cusp][0][0] = a
            matrices[which_cusp][0][1] = b
            matrices[which_cusp][1][0] = c
            matrices[which_cusp][1][1] = d
            if self.is_orientable():
                result = change_peripheral_curves(self.c_triangulation, matrices)
            else:
                result = change_peripheral_curves_nonorientable(
                    self.c_triangulation, matrices)
            free(matrices)
            if result == func_bad_input:
                raise ValueError('The peripheral data ((%d, %d), (%d,%d)) '
                                 'is not acceptable.' % (a,b,c,d))

        else:
            if self.num_cusps() == 1 and len(peripheral_data) == 2:
                self.set_peripheral_curves(peripheral_data, 0)
                return
            if len(peripheral_data) > self.num_cusps():
                raise IndexError('You provided more peripheral data '
                                 'than there are cusps.')

            matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                             sizeof(MatrixInt22))

            for n in range(self.num_cusps()):
                for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                    matrices[n][i][j] = 1 if i == j else 0

            for i, basis in enumerate(peripheral_data):

                meridian, longitude = basis
                a, b = meridian
                c, d = longitude
                if a*d - b*c != 1 and (self.is_orientable() or a*d - b*c != -1):
                    if self.is_orientable():
                        free(matrices)
                        raise ValueError('The data provided does not give a '
                                         '(positively oriented) basis.')
                    else:
                        free(matrices)
                        raise ValueError('The data provided does not give a '
                                         'basis.')

                if 'Klein' in self.cusp_info(i)['topology']:
                    if b != 0 or c != 0:
                        free(matrices)
                        raise ValueError('The data provided are invalid for '
                                         'Klein bottle cusps.')

                matrices[i][0][0] = a
                matrices[i][0][1] = b
                matrices[i][1][0] = c
                matrices[i][1][1] = d

            if self.is_orientable():
                result = change_peripheral_curves(self.c_triangulation, matrices)
            else:
                result = change_peripheral_curves_nonorientable(
                    self.c_triangulation, matrices)
            free(matrices)
            if result == func_bad_input:
                raise ValueError('The peripheral data %s is not acceptable.' %
                                 peripheral_data)

    def has_finite_vertices(self) -> bool:
        """
        Returns ``True`` if and only if the triangulation has finite (non-ideal)
        vertices.

        >>> T = Triangulation("m004")
        >>> T.has_finite_vertices()
        False
        >>> T.dehn_fill((12,13))
        >>> S = T.filled_triangulation()
        >>> S.has_finite_vertices()
        True

        When trying to find a hyperbolic structure, SnapPea will eliminate
        finite vertices:

        >>> M = S.with_hyperbolic_structure()
        >>> M.has_finite_vertices()
        False
        """
        # Bail if empty
        if self.c_triangulation is NULL:
            return False

        count_cusps(self.c_triangulation)
        return get_num_fake_cusps(self.c_triangulation) > 0

    def triangulation_isosig(self,
                             decorated : bool = True,
                             ignore_cusp_ordering : bool = False,
                             ignore_curves : bool = False,
                             ignore_curve_orientations : bool = False,
                             ignore_filling_orientations : bool = False,
                             ignore_orientation : bool = True) -> str:
        """
        Returns the "(decorated) isomorphism signature", a compact text
        representation of the triangulation::

          >>> T = Triangulation('m004')
          >>> T.triangulation_isosig()
          'cPcbbbiht_BaCB'

        This string can be used later to recreate an isomorphic triangulation::

          >>> U = Triangulation('cPcbbbiht_BaCB')
          >>> T == U
          True

        The isomorphism signature is also used to compute the
        :meth:`isometry_signature <snappy.Manifold.isometry_signature>`.
        It comes in two flavors controlled by the :attr:`decorated` flag.

        **Undecorated isomorphism signature**

        The undecorated isomorphism signature is a complete invariant of the
        (oriented) triangulation up to combinatorial isomorphism::

            >>> T = Triangulation('m015')
            >>> T.triangulation_isosig(decorated=False)
            'dLQbcccdero'

        It was introduced in
        `Burton '11 <http://arxiv.org/abs/1110.6080>`_. It canonizes and
        generalizes the ealier dehydration string by
        `Callahan, Hildebrand and Weeks '99
        <https://doi.org/10.1090/S0025-5718-99-01036-4>`_.
        The undecorated isomorphism signature can also be given to
        `Regina <https://regina-normal.github.io/>`_'s
        ``Triangulation3.fromIsoSig``.

        By default, the orientation (if orientable) is ignored. More
        precisely, it computes the string for both orientations (if orientable)
        and uses the lexicographically smaller string::

            >>> T = Triangulation('m015')
            >>> T.triangulation_isosig(decorated=False)
            'dLQbcccdero'
            >>> T.reverse_orientation()
            >>> T.triangulation_isosig(decorated=False)
            'dLQbcccdero'

        When specifying :attr:`ignore_orientation = False`, the result
        encodes the orientation (if orientable). Now the result is
        different if we change the orientation of a chiral triangulation::

            >>> T = Triangulation('m015')
            >>> T.triangulation_isosig(decorated=False, ignore_orientation=False)
            'dLQbcccdero'
            >>> T.reverse_orientation()
            >>> T.triangulation_isosig(decorated=False, ignore_orientation=False)
            'dLQbccceekg'

        **Decorated isomorphism signature (default)**

        SnapPy can decorate the isomorphism signature to include the following
        peripheral information in a canonical way (that is invariant under
        the action by combinatorial isomorphisms of the triangulation):

        #. Indexing of the cusps (that is, ideal vertices).

           * Included by default.
             Can be suppressed with :attr:`ignore_cusp_ordering = True`.

        #. Peripheral curves (aka meridian and longitude, up to homotopy).

           * Included by default.
             Can be suppressed with :attr:`ignore_curves = True`.
           * By default, the decoration encodes the oriented peripheral curves.
             By specifying :attr:`ignore_curve_orientations = True`, it encodes
             the unoriented peripheral curves instead.

        #. Dehn-fillings and orientation of core curve (if present).

           * By default, the decoration encodes the Dehn-fillings and the
             orientation of the core curve.
             That is, we encode the (unoriented) peripheral curve that is
             used for the Dehn-filling as well as the orientation of the class
             of peripheral curves intersecting that peripheral curve ones
             (this explanation only works if the coefficients are co-prime
             integers). Note that this behavior has changed in Version 3.3.
             By specifying :attr:`ignore_filling_orientations = True`, the
             decoration encodes the Dehn-fillings without any regard to
             orientation.
             That is, it normalizes the Dehn-filling coefficients by picking
             a canonical pair among :math:`(m,l)` and :math:`(-m,-l)`.

        Details of the encoding are explained in the
        `SnapPy source code <https://github.com/3-manifolds/SnapPy/blob/master/python/decorated_isosig.py>`_.

        **Example**

        Let us consider the links :math:`9^2_{34}` and ``L9a21``. Note that we use
        :meth:`canonical_retriangulation <snappy.Manifold.canonical_retriangulation>`
        to make the following examples say something intrinsic about the
        hyperbolic manifold::

          >>> from snappy import Manifold
          >>> M = Manifold('9^2_34').canonical_retriangulation()
          >>> N = Manifold('L9a21').canonical_retriangulation()

        The decorated isosig recovers the entire peripheral information faithfully
        (including orientation, see below)::

          >>> M.triangulation_isosig()
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_baBabbbBbC'
          >>> K = Triangulation('oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_baBabbbBbC')
          >>> M.isomorphisms_to(K)
          [0 -> 0  1 -> 1
          [1 0]   [1 0]
          [0 1]   [0 1]
          Extends to link]

        The two links have isometric complements::

          >>> M.triangulation_isosig(decorated=False)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo'
          >>> N.triangulation_isosig(decorated=False)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo'

        However, the complements have different handedness::

          >>> M.triangulation_isosig(decorated=False,ignore_orientation=False)
          'oLLzLPwzQQccdeghjiiklmnmnnuvuvvavovvffffo'
          >>> N.triangulation_isosig(decorated=False,ignore_orientation=False)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo'

        Also, the cusps/components of the link are indexed differently::

          >>> M.triangulation_isosig(ignore_curves=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_ba'
          >>> N.triangulation_isosig(ignore_curves=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_ab'

        Ignoring the indexing, we also see that the oriented merdians and
        longitudes do not match::

          >>> M.triangulation_isosig(ignore_cusp_ordering=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_bBbCBabb'
          >>> N.triangulation_isosig(ignore_cusp_ordering=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_BbbCbabb'

        However, they are the same links (ignoring indexing and orientation)::

          >>> M.triangulation_isosig(ignore_cusp_ordering=True, ignore_curve_orientations=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_bBBcbabb'
          >>> N.triangulation_isosig(ignore_cusp_ordering=True, ignore_curve_orientations=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo_bBBcbabb'

        Let us create two surgery presentations from the links (note that we
        fill after
        :meth:`canonical_retriangulation <snappy.Manifold.canonical_retriangulation>`
        since it rejects Dehn-fillings)::

          >>> M.dehn_fill((4, 5),0)
          >>> N.dehn_fill((4,-5),1)

        They are equivalent surgery presentations (of the same manifold)::

          >>> M.triangulation_isosig(
          ...         ignore_cusp_ordering=True,
          ...         ignore_curves=True,
          ...         ignore_filling_orientations=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo(0,0)(1,5)'
          >>> N.triangulation_isosig(
          ...         ignore_cusp_ordering=True,
          ...         ignore_curves=True,
          ...         ignore_filling_orientations=True)
          'oLLvzQLLQQccdhifihnlmkmlnnpvuvbvouggbggoo(0,0)(1,5)'

        A geometric spun-triangulation of the three-fold cover of m007(3,1)::

          >>> M = Manifold('gLLMQcbeefefujumapn(2,-1)')
          >>> M.solution_type()
          'all tetrahedra positively oriented'

        Specifying :attr:`ignore_filling_orientations = True` can change the
        orientation of the core curve and a triangulation spinning the other
        way that fails to be geometric::

          >>> isosig = M.triangulation_isosig(
          ...         ignore_orientation=False,
          ...         ignore_curves=True,
          ...         ignore_filling_orientations=True,
          ...         ignore_cusp_ordering=True)
          >>> isosig
          'gLLMQcbeefefujumapn(-2,1)'
          >>> Manifold(isosig).solution_type()
          'contains negatively oriented tetrahedra'

        Preserve the unoriented spun-triangulation structure including the
        direction it is spinning::

          >>> isosig = M.triangulation_isosig(
          ...         ignore_cusp_ordering=True,
          ...         ignore_curves=True)
          >>> isosig
          'gLLMQcbeefefpjaqupw(1,1)'
          >>> Manifold(isosig).solution_type()
          'all tetrahedra positively oriented'

        **Orientation**

        Note that :attr:`ignore_orientation=True` only applies to the undecorated
        part of the isomorphism signature. The decoration can still capture the
        the orientation.
        More, precisely, the result of :meth:`.triangulation_isosig` depends on
        the orientation (if the triangulation is orientable and chiral) if any
        of the following is true:

        #. :attr:`ignore_orientation = False`.

        #. :attr:`decorated = True` and
           :attr:`ignore_curves = False` and
           :attr:`ignore_filling_orientations = False`.

        In these cases, re-constructing a triangulation from the isomorphism
        signature yields a triangulation with the same handedness.

        :param decorated:
                Include peripheral information such as indexing of the cusps,
                (oriented or unoriented) peripheral curves and
                (oriented or unoriented) Dehn-fillings.
        :param ignore_cusp_ordering:
                Do not encode the indexing of the cusps.
                Only relevant if :attr:`decorated = True`.
        :param ignore_curves:
                Do not encode the peripheral curves.
                Only relevant if :attr:`decorated = True`.
                This is new in SnapPy version 3.2.
                If :attr:`ignore_curves = True`, the result of this method cannot
                be given to prior versions.
        :param ignore_curve_orientations:
                Do not encode the orientations of the peripheral curves.
                Only relevant if :attr:`decorated = True` and
                :attr:`ignore_curves = False`.
        :param ignore_filling_orientations:
                Do not encode the orientation of the class of peripheral curves
                intersecting the Dehn-fillings.
                Only relevant if :attr:`decorated = True`.
        :param ignore_orientation:
                Do not encode the orientation of the triangulation in the
                undecorated part of the triangulation isosig.
                See above section about orientation.
        """

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        args = (decorated,
                ignore_cusp_ordering,
                ignore_curves,
                ignore_curve_orientations,
                ignore_filling_orientations,
                ignore_orientation)
        try:
            return self._cache.lookup('triangulation_isosig', *args)
        except KeyError:
            pass

        if decorated:
            result = decorated_isosig.decorated_isosig(
                self, _triangulation_class,
                ignore_cusp_ordering = ignore_cusp_ordering,
                ignore_curves = ignore_curves,
                ignore_curve_orientations = ignore_curve_orientations,
                ignore_filling_orientations = ignore_filling_orientations,
                ignore_orientation = ignore_orientation)
        else:
            result = self._undecorated_triangulation_isosig(
                ignore_orientation = ignore_orientation)
        return self._cache.save(result, 'triangulation_isosig', *args)

    def _undecorated_triangulation_isosig(self,
                                          ignore_orientation : bool = False) -> str:
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        cdef char *c_string
        try:
            c_string = get_isomorphism_signature(
                self.c_triangulation, ignore_orientation)
            return to_str(c_string)
        finally:
            free(c_string)

    def _symplectic_form(self, u, v):
        return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i] for i in range(len(u) // 2)])

    def _is_symplectic(self, M):
        """
        Test if the matrix M is symplectic
        :param M: square matrix
        :return: true or false
        """
        n = len(M)

        for i in range(n):
            for j in range(i, n):
                omega = abs(self._symplectic_form(M[i], M[j]))

                if i % 2 == 0 and j % 2 == 1 and j == i + 1:
                    if omega != 2:
                        return False
                elif omega:
                    return False

        return True

    def symplectic_basis(self, verify : bool = False):
        """
        Extend the Neumann-Zagier matrix to one which is symplectic
        (up to factors of 2) using oscillating curves, see
        `Mathews and Purcell '22 <https://arxiv.org/abs/2208.06969>`_.
        Only accepts triangulations with 1 cusp.

          >>> M = Manifold("4_1")
          >>> M.symplectic_basis()
          [-1  0 -1 -1]
          [ 2  0 -2  0]
          [-2 -1 -2 -1]
          [ 0 -1 -2 -1]


        :param verify: Explicitly test if the resulting matrix is symplectic.
        """

        cdef int **c_eqns;
        cdef int **g_eqns;
        cdef int num_rows, num_cols, dual_rows;
        cdef int* eqn;

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        # current get_symplectic_eqns() implementation requires 1 cusp
        if self.num_cusps() > 1:
            raise ValueError('Triangulation contains {} cusps, only accepts triangulations with 1 cusp'.format(self.num_cusps()))

        eqns = []

        peripheral_curves(self.c_triangulation)

        # Cusp Equations
        for i in range(self.num_cusps()):
            cusp_info = self.cusp_info(i)
            if cusp_info.is_complete:
                to_do = [(1,0), (0,1)]
            else:
                to_do = [cusp_info.filling]
            for (m, l) in to_do:
                eqn = get_cusp_equation(self.c_triangulation,
                                        i, int(m), int(l), &num_rows)
                eqns.append([eqn[j] for j in range(num_rows)])
                free_cusp_equation(eqn)

        # Dual Curve Equations
        g_eqns = get_symplectic_basis(self.c_triangulation, &dual_rows, &num_cols, 0)

        for i in range(dual_rows):
            eqns.append([g_eqns[i][j] for j in range(num_cols)])

        free_symplectic_basis(g_eqns, dual_rows)

        # Convert to Neumann Zagier Matrix
        rows = len(eqns)
        retval = [[eqns[i][3 * (j // 2) + j % 2] - eqns[i][3 * (j // 2) + 2] for j in range(rows)] for i in range(rows)]

        if verify:
            if self._is_symplectic(retval):
                print("Result is symplectic (up to factors of 2)")
            else:
                print("Warning: Result is not symplectic")

        return matrix(retval)
