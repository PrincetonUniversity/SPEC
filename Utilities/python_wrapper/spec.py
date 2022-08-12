from __future__ import print_function, absolute_import, division
import _spec
import f90wrap.runtime
import logging

class Inputlist(f90wrap.runtime.FortranModule):
    """
    Module inputlist
    
    
    Defined at inputlist.fpp lines 11-975
    
    """
    @staticmethod
    def initialize_inputs():
        """
        initialize_inputs()
        
        
        Defined at inputlist.fpp lines 845-974
        
        
        """
        _spec.f90wrap_initialize_inputs()
    
    @property
    def mnvol(self):
        """
        Element mnvol ftype=integer pytype=int
        
        
        Defined at inputlist.fpp line 18
        
        """
        return _spec.f90wrap_inputlist__get__mnvol()
    
    @property
    def mmpol(self):
        """
        Element mmpol ftype=integer pytype=int
        
        
        Defined at inputlist.fpp line 19
        
        """
        return _spec.f90wrap_inputlist__get__mmpol()
    
    @property
    def mntor(self):
        """
        Element mntor ftype=integer pytype=int
        
        
        Defined at inputlist.fpp line 20
        
        """
        return _spec.f90wrap_inputlist__get__mntor()
    
    @property
    def igeometry(self):
        """
        Element igeometry ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 23
        
        """
        return _spec.f90wrap_inputlist__get__igeometry()
    
    @igeometry.setter
    def igeometry(self, igeometry):
        _spec.f90wrap_inputlist__set__igeometry(igeometry)
    
    @property
    def istellsym(self):
        """
        Element istellsym ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 24
        
        """
        return _spec.f90wrap_inputlist__get__istellsym()
    
    @istellsym.setter
    def istellsym(self, istellsym):
        _spec.f90wrap_inputlist__set__istellsym(istellsym)
    
    @property
    def lfreebound(self):
        """
        Element lfreebound ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 25
        
        """
        return _spec.f90wrap_inputlist__get__lfreebound()
    
    @lfreebound.setter
    def lfreebound(self, lfreebound):
        _spec.f90wrap_inputlist__set__lfreebound(lfreebound)
    
    @property
    def phiedge(self):
        """
        Element phiedge ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 26
        
        """
        return _spec.f90wrap_inputlist__get__phiedge()
    
    @phiedge.setter
    def phiedge(self, phiedge):
        _spec.f90wrap_inputlist__set__phiedge(phiedge)
    
    @property
    def curtor(self):
        """
        Element curtor ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 27
        
        """
        return _spec.f90wrap_inputlist__get__curtor()
    
    @curtor.setter
    def curtor(self, curtor):
        _spec.f90wrap_inputlist__set__curtor(curtor)
    
    @property
    def curpol(self):
        """
        Element curpol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 28
        
        """
        return _spec.f90wrap_inputlist__get__curpol()
    
    @curpol.setter
    def curpol(self, curpol):
        _spec.f90wrap_inputlist__set__curpol(curpol)
    
    @property
    def gamma(self):
        """
        Element gamma ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 29
        
        """
        return _spec.f90wrap_inputlist__get__gamma()
    
    @gamma.setter
    def gamma(self, gamma):
        _spec.f90wrap_inputlist__set__gamma(gamma)
    
    @property
    def nfp(self):
        """
        Element nfp ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 30
        
        """
        return _spec.f90wrap_inputlist__get__nfp()
    
    @nfp.setter
    def nfp(self, nfp):
        _spec.f90wrap_inputlist__set__nfp(nfp)
    
    @property
    def nvol(self):
        """
        Element nvol ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 31
        
        """
        return _spec.f90wrap_inputlist__get__nvol()
    
    @nvol.setter
    def nvol(self, nvol):
        _spec.f90wrap_inputlist__set__nvol(nvol)
    
    @property
    def mpol(self):
        """
        Element mpol ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 32
        
        """
        return _spec.f90wrap_inputlist__get__mpol()
    
    @mpol.setter
    def mpol(self, mpol):
        _spec.f90wrap_inputlist__set__mpol(mpol)
    
    @property
    def ntor(self):
        """
        Element ntor ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 33
        
        """
        return _spec.f90wrap_inputlist__get__ntor()
    
    @ntor.setter
    def ntor(self, ntor):
        _spec.f90wrap_inputlist__set__ntor(ntor)
    
    @property
    def lrad(self):
        """
        Element lrad ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__lrad(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lrad = self._arrays[array_handle]
        else:
            lrad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__lrad)
            self._arrays[array_handle] = lrad
        return lrad
    
    @lrad.setter
    def lrad(self, lrad):
        self.lrad[...] = lrad
    
    @property
    def lconstraint(self):
        """
        Element lconstraint ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 35
        
        """
        return _spec.f90wrap_inputlist__get__lconstraint()
    
    @lconstraint.setter
    def lconstraint(self, lconstraint):
        _spec.f90wrap_inputlist__set__lconstraint(lconstraint)
    
    @property
    def tflux(self):
        """
        Element tflux ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__tflux(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tflux = self._arrays[array_handle]
        else:
            tflux = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__tflux)
            self._arrays[array_handle] = tflux
        return tflux
    
    @tflux.setter
    def tflux(self, tflux):
        self.tflux[...] = tflux
    
    @property
    def pflux(self):
        """
        Element pflux ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__pflux(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pflux = self._arrays[array_handle]
        else:
            pflux = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__pflux)
            self._arrays[array_handle] = pflux
        return pflux
    
    @pflux.setter
    def pflux(self, pflux):
        self.pflux[...] = pflux
    
    @property
    def helicity(self):
        """
        Element helicity ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__helicity(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            helicity = self._arrays[array_handle]
        else:
            helicity = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__helicity)
            self._arrays[array_handle] = helicity
        return helicity
    
    @helicity.setter
    def helicity(self, helicity):
        self.helicity[...] = helicity
    
    @property
    def pscale(self):
        """
        Element pscale ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 39
        
        """
        return _spec.f90wrap_inputlist__get__pscale()
    
    @pscale.setter
    def pscale(self, pscale):
        _spec.f90wrap_inputlist__set__pscale(pscale)
    
    @property
    def pressure(self):
        """
        Element pressure ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__pressure(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pressure = self._arrays[array_handle]
        else:
            pressure = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__pressure)
            self._arrays[array_handle] = pressure
        return pressure
    
    @pressure.setter
    def pressure(self, pressure):
        self.pressure[...] = pressure
    
    @property
    def ladiabatic(self):
        """
        Element ladiabatic ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 41
        
        """
        return _spec.f90wrap_inputlist__get__ladiabatic()
    
    @ladiabatic.setter
    def ladiabatic(self, ladiabatic):
        _spec.f90wrap_inputlist__set__ladiabatic(ladiabatic)
    
    @property
    def adiabatic(self):
        """
        Element adiabatic ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__adiabatic(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            adiabatic = self._arrays[array_handle]
        else:
            adiabatic = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__adiabatic)
            self._arrays[array_handle] = adiabatic
        return adiabatic
    
    @adiabatic.setter
    def adiabatic(self, adiabatic):
        self.adiabatic[...] = adiabatic
    
    @property
    def mu(self):
        """
        Element mu ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__mu(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            mu = self._arrays[array_handle]
        else:
            mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__mu)
            self._arrays[array_handle] = mu
        return mu
    
    @mu.setter
    def mu(self, mu):
        self.mu[...] = mu
    
    @property
    def ivolume(self):
        """
        Element ivolume ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__ivolume(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ivolume = self._arrays[array_handle]
        else:
            ivolume = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__ivolume)
            self._arrays[array_handle] = ivolume
        return ivolume
    
    @ivolume.setter
    def ivolume(self, ivolume):
        self.ivolume[...] = ivolume
    
    @property
    def isurf(self):
        """
        Element isurf ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 45
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__isurf(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            isurf = self._arrays[array_handle]
        else:
            isurf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__isurf)
            self._arrays[array_handle] = isurf
        return isurf
    
    @isurf.setter
    def isurf(self, isurf):
        self.isurf[...] = isurf
    
    @property
    def pl(self):
        """
        Element pl ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__pl(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pl = self._arrays[array_handle]
        else:
            pl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__pl)
            self._arrays[array_handle] = pl
        return pl
    
    @pl.setter
    def pl(self, pl):
        self.pl[...] = pl
    
    @property
    def ql(self):
        """
        Element ql ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 47
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__ql(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ql = self._arrays[array_handle]
        else:
            ql = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__ql)
            self._arrays[array_handle] = ql
        return ql
    
    @ql.setter
    def ql(self, ql):
        self.ql[...] = ql
    
    @property
    def pr(self):
        """
        Element pr ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 48
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__pr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pr = self._arrays[array_handle]
        else:
            pr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__pr)
            self._arrays[array_handle] = pr
        return pr
    
    @pr.setter
    def pr(self, pr):
        self.pr[...] = pr
    
    @property
    def qr(self):
        """
        Element qr ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__qr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            qr = self._arrays[array_handle]
        else:
            qr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__qr)
            self._arrays[array_handle] = qr
        return qr
    
    @qr.setter
    def qr(self, qr):
        self.qr[...] = qr
    
    @property
    def iota(self):
        """
        Element iota ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 50
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__iota(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iota = self._arrays[array_handle]
        else:
            iota = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__iota)
            self._arrays[array_handle] = iota
        return iota
    
    @iota.setter
    def iota(self, iota):
        self.iota[...] = iota
    
    @property
    def lp(self):
        """
        Element lp ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__lp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lp = self._arrays[array_handle]
        else:
            lp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__lp)
            self._arrays[array_handle] = lp
        return lp
    
    @lp.setter
    def lp(self, lp):
        self.lp[...] = lp
    
    @property
    def lq(self):
        """
        Element lq ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__lq(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lq = self._arrays[array_handle]
        else:
            lq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__lq)
            self._arrays[array_handle] = lq
        return lq
    
    @lq.setter
    def lq(self, lq):
        self.lq[...] = lq
    
    @property
    def rp(self):
        """
        Element rp ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rp = self._arrays[array_handle]
        else:
            rp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rp)
            self._arrays[array_handle] = rp
        return rp
    
    @rp.setter
    def rp(self, rp):
        self.rp[...] = rp
    
    @property
    def rq(self):
        """
        Element rq ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rq(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rq = self._arrays[array_handle]
        else:
            rq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rq)
            self._arrays[array_handle] = rq
        return rq
    
    @rq.setter
    def rq(self, rq):
        self.rq[...] = rq
    
    @property
    def oita(self):
        """
        Element oita ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 55
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__oita(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            oita = self._arrays[array_handle]
        else:
            oita = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__oita)
            self._arrays[array_handle] = oita
        return oita
    
    @oita.setter
    def oita(self, oita):
        self.oita[...] = oita
    
    @property
    def rpol(self):
        """
        Element rpol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 56
        
        """
        return _spec.f90wrap_inputlist__get__rpol()
    
    @rpol.setter
    def rpol(self, rpol):
        _spec.f90wrap_inputlist__set__rpol(rpol)
    
    @property
    def rtor(self):
        """
        Element rtor ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 57
        
        """
        return _spec.f90wrap_inputlist__get__rtor()
    
    @rtor.setter
    def rtor(self, rtor):
        _spec.f90wrap_inputlist__set__rtor(rtor)
    
    @property
    def rac(self):
        """
        Element rac ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 58
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rac(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rac = self._arrays[array_handle]
        else:
            rac = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rac)
            self._arrays[array_handle] = rac
        return rac
    
    @rac.setter
    def rac(self, rac):
        self.rac[...] = rac
    
    @property
    def zas(self):
        """
        Element zas ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 59
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__zas(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zas = self._arrays[array_handle]
        else:
            zas = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__zas)
            self._arrays[array_handle] = zas
        return zas
    
    @zas.setter
    def zas(self, zas):
        self.zas[...] = zas
    
    @property
    def ras(self):
        """
        Element ras ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 60
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__ras(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ras = self._arrays[array_handle]
        else:
            ras = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__ras)
            self._arrays[array_handle] = ras
        return ras
    
    @ras.setter
    def ras(self, ras):
        self.ras[...] = ras
    
    @property
    def zac(self):
        """
        Element zac ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 61
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__zac(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zac = self._arrays[array_handle]
        else:
            zac = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__zac)
            self._arrays[array_handle] = zac
        return zac
    
    @zac.setter
    def zac(self, zac):
        self.zac[...] = zac
    
    @property
    def rbc(self):
        """
        Element rbc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 62
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rbc = self._arrays[array_handle]
        else:
            rbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rbc)
            self._arrays[array_handle] = rbc
        return rbc
    
    @rbc.setter
    def rbc(self, rbc):
        self.rbc[...] = rbc
    
    @property
    def zbs(self):
        """
        Element zbs ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__zbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zbs = self._arrays[array_handle]
        else:
            zbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__zbs)
            self._arrays[array_handle] = zbs
        return zbs
    
    @zbs.setter
    def zbs(self, zbs):
        self.zbs[...] = zbs
    
    @property
    def rbs(self):
        """
        Element rbs ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rbs = self._arrays[array_handle]
        else:
            rbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rbs)
            self._arrays[array_handle] = rbs
        return rbs
    
    @rbs.setter
    def rbs(self, rbs):
        self.rbs[...] = rbs
    
    @property
    def zbc(self):
        """
        Element zbc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 65
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__zbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zbc = self._arrays[array_handle]
        else:
            zbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__zbc)
            self._arrays[array_handle] = zbc
        return zbc
    
    @zbc.setter
    def zbc(self, zbc):
        self.zbc[...] = zbc
    
    @property
    def rwc(self):
        """
        Element rwc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 66
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rwc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rwc = self._arrays[array_handle]
        else:
            rwc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rwc)
            self._arrays[array_handle] = rwc
        return rwc
    
    @rwc.setter
    def rwc(self, rwc):
        self.rwc[...] = rwc
    
    @property
    def zws(self):
        """
        Element zws ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 67
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__zws(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zws = self._arrays[array_handle]
        else:
            zws = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__zws)
            self._arrays[array_handle] = zws
        return zws
    
    @zws.setter
    def zws(self, zws):
        self.zws[...] = zws
    
    @property
    def rws(self):
        """
        Element rws ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 68
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__rws(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rws = self._arrays[array_handle]
        else:
            rws = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__rws)
            self._arrays[array_handle] = rws
        return rws
    
    @rws.setter
    def rws(self, rws):
        self.rws[...] = rws
    
    @property
    def zwc(self):
        """
        Element zwc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 69
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__zwc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zwc = self._arrays[array_handle]
        else:
            zwc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__zwc)
            self._arrays[array_handle] = zwc
        return zwc
    
    @zwc.setter
    def zwc(self, zwc):
        self.zwc[...] = zwc
    
    @property
    def vns(self):
        """
        Element vns ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 70
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__vns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            vns = self._arrays[array_handle]
        else:
            vns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__vns)
            self._arrays[array_handle] = vns
        return vns
    
    @vns.setter
    def vns(self, vns):
        self.vns[...] = vns
    
    @property
    def bns(self):
        """
        Element bns ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 71
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__bns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bns = self._arrays[array_handle]
        else:
            bns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__bns)
            self._arrays[array_handle] = bns
        return bns
    
    @bns.setter
    def bns(self, bns):
        self.bns[...] = bns
    
    @property
    def vnc(self):
        """
        Element vnc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 72
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__vnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            vnc = self._arrays[array_handle]
        else:
            vnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__vnc)
            self._arrays[array_handle] = vnc
        return vnc
    
    @vnc.setter
    def vnc(self, vnc):
        self.vnc[...] = vnc
    
    @property
    def bnc(self):
        """
        Element bnc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 73
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__bnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bnc = self._arrays[array_handle]
        else:
            bnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__bnc)
            self._arrays[array_handle] = bnc
        return bnc
    
    @bnc.setter
    def bnc(self, bnc):
        self.bnc[...] = bnc
    
    @property
    def mupftol(self):
        """
        Element mupftol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 74
        
        """
        return _spec.f90wrap_inputlist__get__mupftol()
    
    @mupftol.setter
    def mupftol(self, mupftol):
        _spec.f90wrap_inputlist__set__mupftol(mupftol)
    
    @property
    def mupfits(self):
        """
        Element mupfits ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 75
        
        """
        return _spec.f90wrap_inputlist__get__mupfits()
    
    @mupfits.setter
    def mupfits(self, mupfits):
        _spec.f90wrap_inputlist__set__mupfits(mupfits)
    
    @property
    def lreflect(self):
        """
        Element lreflect ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 76
        
        """
        return _spec.f90wrap_inputlist__get__lreflect()
    
    @lreflect.setter
    def lreflect(self, lreflect):
        _spec.f90wrap_inputlist__set__lreflect(lreflect)
    
    @property
    def linitialize(self):
        """
        Element linitialize ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 79
        
        """
        return _spec.f90wrap_inputlist__get__linitialize()
    
    @linitialize.setter
    def linitialize(self, linitialize):
        _spec.f90wrap_inputlist__set__linitialize(linitialize)
    
    @property
    def lautoinitbn(self):
        """
        Element lautoinitbn ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 80
        
        """
        return _spec.f90wrap_inputlist__get__lautoinitbn()
    
    @lautoinitbn.setter
    def lautoinitbn(self, lautoinitbn):
        _spec.f90wrap_inputlist__set__lautoinitbn(lautoinitbn)
    
    @property
    def lzerovac(self):
        """
        Element lzerovac ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 81
        
        """
        return _spec.f90wrap_inputlist__get__lzerovac()
    
    @lzerovac.setter
    def lzerovac(self, lzerovac):
        _spec.f90wrap_inputlist__set__lzerovac(lzerovac)
    
    @property
    def ndiscrete(self):
        """
        Element ndiscrete ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 82
        
        """
        return _spec.f90wrap_inputlist__get__ndiscrete()
    
    @ndiscrete.setter
    def ndiscrete(self, ndiscrete):
        _spec.f90wrap_inputlist__set__ndiscrete(ndiscrete)
    
    @property
    def nquad(self):
        """
        Element nquad ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 83
        
        """
        return _spec.f90wrap_inputlist__get__nquad()
    
    @nquad.setter
    def nquad(self, nquad):
        _spec.f90wrap_inputlist__set__nquad(nquad)
    
    @property
    def impol(self):
        """
        Element impol ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 84
        
        """
        return _spec.f90wrap_inputlist__get__impol()
    
    @impol.setter
    def impol(self, impol):
        _spec.f90wrap_inputlist__set__impol(impol)
    
    @property
    def intor(self):
        """
        Element intor ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 85
        
        """
        return _spec.f90wrap_inputlist__get__intor()
    
    @intor.setter
    def intor(self, intor):
        _spec.f90wrap_inputlist__set__intor(intor)
    
    @property
    def lsparse(self):
        """
        Element lsparse ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 86
        
        """
        return _spec.f90wrap_inputlist__get__lsparse()
    
    @lsparse.setter
    def lsparse(self, lsparse):
        _spec.f90wrap_inputlist__set__lsparse(lsparse)
    
    @property
    def lsvdiota(self):
        """
        Element lsvdiota ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 87
        
        """
        return _spec.f90wrap_inputlist__get__lsvdiota()
    
    @lsvdiota.setter
    def lsvdiota(self, lsvdiota):
        _spec.f90wrap_inputlist__set__lsvdiota(lsvdiota)
    
    @property
    def imethod(self):
        """
        Element imethod ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 88
        
        """
        return _spec.f90wrap_inputlist__get__imethod()
    
    @imethod.setter
    def imethod(self, imethod):
        _spec.f90wrap_inputlist__set__imethod(imethod)
    
    @property
    def iorder(self):
        """
        Element iorder ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 89
        
        """
        return _spec.f90wrap_inputlist__get__iorder()
    
    @iorder.setter
    def iorder(self, iorder):
        _spec.f90wrap_inputlist__set__iorder(iorder)
    
    @property
    def iprecon(self):
        """
        Element iprecon ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 90
        
        """
        return _spec.f90wrap_inputlist__get__iprecon()
    
    @iprecon.setter
    def iprecon(self, iprecon):
        _spec.f90wrap_inputlist__set__iprecon(iprecon)
    
    @property
    def iotatol(self):
        """
        Element iotatol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 91
        
        """
        return _spec.f90wrap_inputlist__get__iotatol()
    
    @iotatol.setter
    def iotatol(self, iotatol):
        _spec.f90wrap_inputlist__set__iotatol(iotatol)
    
    @property
    def lextrap(self):
        """
        Element lextrap ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 92
        
        """
        return _spec.f90wrap_inputlist__get__lextrap()
    
    @lextrap.setter
    def lextrap(self, lextrap):
        _spec.f90wrap_inputlist__set__lextrap(lextrap)
    
    @property
    def mregular(self):
        """
        Element mregular ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 93
        
        """
        return _spec.f90wrap_inputlist__get__mregular()
    
    @mregular.setter
    def mregular(self, mregular):
        _spec.f90wrap_inputlist__set__mregular(mregular)
    
    @property
    def lrzaxis(self):
        """
        Element lrzaxis ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 94
        
        """
        return _spec.f90wrap_inputlist__get__lrzaxis()
    
    @lrzaxis.setter
    def lrzaxis(self, lrzaxis):
        _spec.f90wrap_inputlist__set__lrzaxis(lrzaxis)
    
    @property
    def ntoraxis(self):
        """
        Element ntoraxis ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 95
        
        """
        return _spec.f90wrap_inputlist__get__ntoraxis()
    
    @ntoraxis.setter
    def ntoraxis(self, ntoraxis):
        _spec.f90wrap_inputlist__set__ntoraxis(ntoraxis)
    
    @property
    def lbeltrami(self):
        """
        Element lbeltrami ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 98
        
        """
        return _spec.f90wrap_inputlist__get__lbeltrami()
    
    @lbeltrami.setter
    def lbeltrami(self, lbeltrami):
        _spec.f90wrap_inputlist__set__lbeltrami(lbeltrami)
    
    @property
    def linitgues(self):
        """
        Element linitgues ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 99
        
        """
        return _spec.f90wrap_inputlist__get__linitgues()
    
    @linitgues.setter
    def linitgues(self, linitgues):
        _spec.f90wrap_inputlist__set__linitgues(linitgues)
    
    @property
    def lposdef(self):
        """
        Element lposdef ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 100
        
        """
        return _spec.f90wrap_inputlist__get__lposdef()
    
    @lposdef.setter
    def lposdef(self, lposdef):
        _spec.f90wrap_inputlist__set__lposdef(lposdef)
    
    @property
    def maxrndgues(self):
        """
        Element maxrndgues ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 101
        
        """
        return _spec.f90wrap_inputlist__get__maxrndgues()
    
    @maxrndgues.setter
    def maxrndgues(self, maxrndgues):
        _spec.f90wrap_inputlist__set__maxrndgues(maxrndgues)
    
    @property
    def lmatsolver(self):
        """
        Element lmatsolver ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 102
        
        """
        return _spec.f90wrap_inputlist__get__lmatsolver()
    
    @lmatsolver.setter
    def lmatsolver(self, lmatsolver):
        _spec.f90wrap_inputlist__set__lmatsolver(lmatsolver)
    
    @property
    def nitergmres(self):
        """
        Element nitergmres ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 103
        
        """
        return _spec.f90wrap_inputlist__get__nitergmres()
    
    @nitergmres.setter
    def nitergmres(self, nitergmres):
        _spec.f90wrap_inputlist__set__nitergmres(nitergmres)
    
    @property
    def epsgmres(self):
        """
        Element epsgmres ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 104
        
        """
        return _spec.f90wrap_inputlist__get__epsgmres()
    
    @epsgmres.setter
    def epsgmres(self, epsgmres):
        _spec.f90wrap_inputlist__set__epsgmres(epsgmres)
    
    @property
    def lgmresprec(self):
        """
        Element lgmresprec ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 105
        
        """
        return _spec.f90wrap_inputlist__get__lgmresprec()
    
    @lgmresprec.setter
    def lgmresprec(self, lgmresprec):
        _spec.f90wrap_inputlist__set__lgmresprec(lgmresprec)
    
    @property
    def epsilu(self):
        """
        Element epsilu ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 106
        
        """
        return _spec.f90wrap_inputlist__get__epsilu()
    
    @epsilu.setter
    def epsilu(self, epsilu):
        _spec.f90wrap_inputlist__set__epsilu(epsilu)
    
    @property
    def lfindzero(self):
        """
        Element lfindzero ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 109
        
        """
        return _spec.f90wrap_inputlist__get__lfindzero()
    
    @lfindzero.setter
    def lfindzero(self, lfindzero):
        _spec.f90wrap_inputlist__set__lfindzero(lfindzero)
    
    @property
    def escale(self):
        """
        Element escale ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 110
        
        """
        return _spec.f90wrap_inputlist__get__escale()
    
    @escale.setter
    def escale(self, escale):
        _spec.f90wrap_inputlist__set__escale(escale)
    
    @property
    def opsilon(self):
        """
        Element opsilon ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 111
        
        """
        return _spec.f90wrap_inputlist__get__opsilon()
    
    @opsilon.setter
    def opsilon(self, opsilon):
        _spec.f90wrap_inputlist__set__opsilon(opsilon)
    
    @property
    def pcondense(self):
        """
        Element pcondense ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 112
        
        """
        return _spec.f90wrap_inputlist__get__pcondense()
    
    @pcondense.setter
    def pcondense(self, pcondense):
        _spec.f90wrap_inputlist__set__pcondense(pcondense)
    
    @property
    def epsilon(self):
        """
        Element epsilon ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 113
        
        """
        return _spec.f90wrap_inputlist__get__epsilon()
    
    @epsilon.setter
    def epsilon(self, epsilon):
        _spec.f90wrap_inputlist__set__epsilon(epsilon)
    
    @property
    def wpoloidal(self):
        """
        Element wpoloidal ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 114
        
        """
        return _spec.f90wrap_inputlist__get__wpoloidal()
    
    @wpoloidal.setter
    def wpoloidal(self, wpoloidal):
        _spec.f90wrap_inputlist__set__wpoloidal(wpoloidal)
    
    @property
    def upsilon(self):
        """
        Element upsilon ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 115
        
        """
        return _spec.f90wrap_inputlist__get__upsilon()
    
    @upsilon.setter
    def upsilon(self, upsilon):
        _spec.f90wrap_inputlist__set__upsilon(upsilon)
    
    @property
    def forcetol(self):
        """
        Element forcetol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 116
        
        """
        return _spec.f90wrap_inputlist__get__forcetol()
    
    @forcetol.setter
    def forcetol(self, forcetol):
        _spec.f90wrap_inputlist__set__forcetol(forcetol)
    
    @property
    def c05xmax(self):
        """
        Element c05xmax ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 117
        
        """
        return _spec.f90wrap_inputlist__get__c05xmax()
    
    @c05xmax.setter
    def c05xmax(self, c05xmax):
        _spec.f90wrap_inputlist__set__c05xmax(c05xmax)
    
    @property
    def c05xtol(self):
        """
        Element c05xtol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 118
        
        """
        return _spec.f90wrap_inputlist__get__c05xtol()
    
    @c05xtol.setter
    def c05xtol(self, c05xtol):
        _spec.f90wrap_inputlist__set__c05xtol(c05xtol)
    
    @property
    def c05factor(self):
        """
        Element c05factor ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 119
        
        """
        return _spec.f90wrap_inputlist__get__c05factor()
    
    @c05factor.setter
    def c05factor(self, c05factor):
        _spec.f90wrap_inputlist__set__c05factor(c05factor)
    
    @property
    def lreadgf(self):
        """
        Element lreadgf ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 120
        
        """
        return _spec.f90wrap_inputlist__get__lreadgf()
    
    @lreadgf.setter
    def lreadgf(self, lreadgf):
        _spec.f90wrap_inputlist__set__lreadgf(lreadgf)
    
    @property
    def mfreeits(self):
        """
        Element mfreeits ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 121
        
        """
        return _spec.f90wrap_inputlist__get__mfreeits()
    
    @mfreeits.setter
    def mfreeits(self, mfreeits):
        _spec.f90wrap_inputlist__set__mfreeits(mfreeits)
    
    @property
    def bnstol(self):
        """
        Element bnstol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 122
        
        """
        return _spec.f90wrap_inputlist__get__bnstol()
    
    @bnstol.setter
    def bnstol(self, bnstol):
        _spec.f90wrap_inputlist__set__bnstol(bnstol)
    
    @property
    def bnsblend(self):
        """
        Element bnsblend ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 123
        
        """
        return _spec.f90wrap_inputlist__get__bnsblend()
    
    @bnsblend.setter
    def bnsblend(self, bnsblend):
        _spec.f90wrap_inputlist__set__bnsblend(bnsblend)
    
    @property
    def gbntol(self):
        """
        Element gbntol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 124
        
        """
        return _spec.f90wrap_inputlist__get__gbntol()
    
    @gbntol.setter
    def gbntol(self, gbntol):
        _spec.f90wrap_inputlist__set__gbntol(gbntol)
    
    @property
    def gbnbld(self):
        """
        Element gbnbld ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 125
        
        """
        return _spec.f90wrap_inputlist__get__gbnbld()
    
    @gbnbld.setter
    def gbnbld(self, gbnbld):
        _spec.f90wrap_inputlist__set__gbnbld(gbnbld)
    
    @property
    def vcasingeps(self):
        """
        Element vcasingeps ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 126
        
        """
        return _spec.f90wrap_inputlist__get__vcasingeps()
    
    @vcasingeps.setter
    def vcasingeps(self, vcasingeps):
        _spec.f90wrap_inputlist__set__vcasingeps(vcasingeps)
    
    @property
    def vcasingtol(self):
        """
        Element vcasingtol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 127
        
        """
        return _spec.f90wrap_inputlist__get__vcasingtol()
    
    @vcasingtol.setter
    def vcasingtol(self, vcasingtol):
        _spec.f90wrap_inputlist__set__vcasingtol(vcasingtol)
    
    @property
    def vcasingits(self):
        """
        Element vcasingits ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 128
        
        """
        return _spec.f90wrap_inputlist__get__vcasingits()
    
    @vcasingits.setter
    def vcasingits(self, vcasingits):
        _spec.f90wrap_inputlist__set__vcasingits(vcasingits)
    
    @property
    def vcasingper(self):
        """
        Element vcasingper ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 129
        
        """
        return _spec.f90wrap_inputlist__get__vcasingper()
    
    @vcasingper.setter
    def vcasingper(self, vcasingper):
        _spec.f90wrap_inputlist__set__vcasingper(vcasingper)
    
    @property
    def mcasingcal(self):
        """
        Element mcasingcal ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 130
        
        """
        return _spec.f90wrap_inputlist__get__mcasingcal()
    
    @mcasingcal.setter
    def mcasingcal(self, mcasingcal):
        _spec.f90wrap_inputlist__set__mcasingcal(mcasingcal)
    
    @property
    def odetol(self):
        """
        Element odetol ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 133
        
        """
        return _spec.f90wrap_inputlist__get__odetol()
    
    @odetol.setter
    def odetol(self, odetol):
        _spec.f90wrap_inputlist__set__odetol(odetol)
    
    @property
    def absreq(self):
        """
        Element absreq ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 134
        
        """
        return _spec.f90wrap_inputlist__get__absreq()
    
    @absreq.setter
    def absreq(self, absreq):
        _spec.f90wrap_inputlist__set__absreq(absreq)
    
    @property
    def relreq(self):
        """
        Element relreq ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 135
        
        """
        return _spec.f90wrap_inputlist__get__relreq()
    
    @relreq.setter
    def relreq(self, relreq):
        _spec.f90wrap_inputlist__set__relreq(relreq)
    
    @property
    def absacc(self):
        """
        Element absacc ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 136
        
        """
        return _spec.f90wrap_inputlist__get__absacc()
    
    @absacc.setter
    def absacc(self, absacc):
        _spec.f90wrap_inputlist__set__absacc(absacc)
    
    @property
    def epsr(self):
        """
        Element epsr ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 137
        
        """
        return _spec.f90wrap_inputlist__get__epsr()
    
    @epsr.setter
    def epsr(self, epsr):
        _spec.f90wrap_inputlist__set__epsr(epsr)
    
    @property
    def nppts(self):
        """
        Element nppts ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 138
        
        """
        return _spec.f90wrap_inputlist__get__nppts()
    
    @nppts.setter
    def nppts(self, nppts):
        _spec.f90wrap_inputlist__set__nppts(nppts)
    
    @property
    def ppts(self):
        """
        Element ppts ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 139
        
        """
        return _spec.f90wrap_inputlist__get__ppts()
    
    @ppts.setter
    def ppts(self, ppts):
        _spec.f90wrap_inputlist__set__ppts(ppts)
    
    @property
    def nptrj(self):
        """
        Element nptrj ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 140
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_inputlist__array__nptrj(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            nptrj = self._arrays[array_handle]
        else:
            nptrj = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_inputlist__array__nptrj)
            self._arrays[array_handle] = nptrj
        return nptrj
    
    @nptrj.setter
    def nptrj(self, nptrj):
        self.nptrj[...] = nptrj
    
    @property
    def lhevalues(self):
        """
        Element lhevalues ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 141
        
        """
        return _spec.f90wrap_inputlist__get__lhevalues()
    
    @lhevalues.setter
    def lhevalues(self, lhevalues):
        _spec.f90wrap_inputlist__set__lhevalues(lhevalues)
    
    @property
    def lhevectors(self):
        """
        Element lhevectors ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 142
        
        """
        return _spec.f90wrap_inputlist__get__lhevectors()
    
    @lhevectors.setter
    def lhevectors(self, lhevectors):
        _spec.f90wrap_inputlist__set__lhevectors(lhevectors)
    
    @property
    def lhmatrix(self):
        """
        Element lhmatrix ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 143
        
        """
        return _spec.f90wrap_inputlist__get__lhmatrix()
    
    @lhmatrix.setter
    def lhmatrix(self, lhmatrix):
        _spec.f90wrap_inputlist__set__lhmatrix(lhmatrix)
    
    @property
    def lperturbed(self):
        """
        Element lperturbed ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 144
        
        """
        return _spec.f90wrap_inputlist__get__lperturbed()
    
    @lperturbed.setter
    def lperturbed(self, lperturbed):
        _spec.f90wrap_inputlist__set__lperturbed(lperturbed)
    
    @property
    def dpp(self):
        """
        Element dpp ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 145
        
        """
        return _spec.f90wrap_inputlist__get__dpp()
    
    @dpp.setter
    def dpp(self, dpp):
        _spec.f90wrap_inputlist__set__dpp(dpp)
    
    @property
    def dqq(self):
        """
        Element dqq ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 146
        
        """
        return _spec.f90wrap_inputlist__get__dqq()
    
    @dqq.setter
    def dqq(self, dqq):
        _spec.f90wrap_inputlist__set__dqq(dqq)
    
    @property
    def lerrortype(self):
        """
        Element lerrortype ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 147
        
        """
        return _spec.f90wrap_inputlist__get__lerrortype()
    
    @lerrortype.setter
    def lerrortype(self, lerrortype):
        _spec.f90wrap_inputlist__set__lerrortype(lerrortype)
    
    @property
    def ngrid(self):
        """
        Element ngrid ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 148
        
        """
        return _spec.f90wrap_inputlist__get__ngrid()
    
    @ngrid.setter
    def ngrid(self, ngrid):
        _spec.f90wrap_inputlist__set__ngrid(ngrid)
    
    @property
    def drz(self):
        """
        Element drz ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 149
        
        """
        return _spec.f90wrap_inputlist__get__drz()
    
    @drz.setter
    def drz(self, drz):
        _spec.f90wrap_inputlist__set__drz(drz)
    
    @property
    def lcheck(self):
        """
        Element lcheck ftype=integer       pytype=int
        
        
        Defined at inputlist.fpp line 150
        
        """
        return _spec.f90wrap_inputlist__get__lcheck()
    
    @lcheck.setter
    def lcheck(self, lcheck):
        _spec.f90wrap_inputlist__set__lcheck(lcheck)
    
    @property
    def ltiming(self):
        """
        Element ltiming ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 151
        
        """
        return _spec.f90wrap_inputlist__get__ltiming()
    
    @ltiming.setter
    def ltiming(self, ltiming):
        _spec.f90wrap_inputlist__set__ltiming(ltiming)
    
    @property
    def fudge(self):
        """
        Element fudge ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 152
        
        """
        return _spec.f90wrap_inputlist__get__fudge()
    
    @fudge.setter
    def fudge(self, fudge):
        _spec.f90wrap_inputlist__set__fudge(fudge)
    
    @property
    def scaling(self):
        """
        Element scaling ftype=real(8) pytype=float
        
        
        Defined at inputlist.fpp line 153
        
        """
        return _spec.f90wrap_inputlist__get__scaling()
    
    @scaling.setter
    def scaling(self, scaling):
        _spec.f90wrap_inputlist__set__scaling(scaling)
    
    @property
    def wdcuhre(self):
        """
        Element wdcuhre ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 156
        
        """
        return _spec.f90wrap_inputlist__get__wdcuhre()
    
    @wdcuhre.setter
    def wdcuhre(self, wdcuhre):
        _spec.f90wrap_inputlist__set__wdcuhre(wdcuhre)
    
    @property
    def wminpack(self):
        """
        Element wminpack ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 157
        
        """
        return _spec.f90wrap_inputlist__get__wminpack()
    
    @wminpack.setter
    def wminpack(self, wminpack):
        _spec.f90wrap_inputlist__set__wminpack(wminpack)
    
    @property
    def wiqpack(self):
        """
        Element wiqpack ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 158
        
        """
        return _spec.f90wrap_inputlist__get__wiqpack()
    
    @wiqpack.setter
    def wiqpack(self, wiqpack):
        _spec.f90wrap_inputlist__set__wiqpack(wiqpack)
    
    @property
    def wrksuite(self):
        """
        Element wrksuite ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 159
        
        """
        return _spec.f90wrap_inputlist__get__wrksuite()
    
    @wrksuite.setter
    def wrksuite(self, wrksuite):
        _spec.f90wrap_inputlist__set__wrksuite(wrksuite)
    
    @property
    def wi1mach(self):
        """
        Element wi1mach ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 160
        
        """
        return _spec.f90wrap_inputlist__get__wi1mach()
    
    @wi1mach.setter
    def wi1mach(self, wi1mach):
        _spec.f90wrap_inputlist__set__wi1mach(wi1mach)
    
    @property
    def wd1mach(self):
        """
        Element wd1mach ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 161
        
        """
        return _spec.f90wrap_inputlist__get__wd1mach()
    
    @wd1mach.setter
    def wd1mach(self, wd1mach):
        _spec.f90wrap_inputlist__set__wd1mach(wd1mach)
    
    @property
    def wilut(self):
        """
        Element wilut ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 162
        
        """
        return _spec.f90wrap_inputlist__get__wilut()
    
    @wilut.setter
    def wilut(self, wilut):
        _spec.f90wrap_inputlist__set__wilut(wilut)
    
    @property
    def witers(self):
        """
        Element witers ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 163
        
        """
        return _spec.f90wrap_inputlist__get__witers()
    
    @witers.setter
    def witers(self, witers):
        _spec.f90wrap_inputlist__set__witers(witers)
    
    @property
    def winputlist(self):
        """
        Element winputlist ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 164
        
        """
        return _spec.f90wrap_inputlist__get__winputlist()
    
    @winputlist.setter
    def winputlist(self, winputlist):
        _spec.f90wrap_inputlist__set__winputlist(winputlist)
    
    @property
    def wglobal(self):
        """
        Element wglobal ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 165
        
        """
        return _spec.f90wrap_inputlist__get__wglobal()
    
    @wglobal.setter
    def wglobal(self, wglobal):
        _spec.f90wrap_inputlist__set__wglobal(wglobal)
    
    @property
    def wsphdf5(self):
        """
        Element wsphdf5 ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 166
        
        """
        return _spec.f90wrap_inputlist__get__wsphdf5()
    
    @wsphdf5.setter
    def wsphdf5(self, wsphdf5):
        _spec.f90wrap_inputlist__set__wsphdf5(wsphdf5)
    
    @property
    def wpreset(self):
        """
        Element wpreset ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 167
        
        """
        return _spec.f90wrap_inputlist__get__wpreset()
    
    @wpreset.setter
    def wpreset(self, wpreset):
        _spec.f90wrap_inputlist__set__wpreset(wpreset)
    
    @property
    def wmanual(self):
        """
        Element wmanual ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 168
        
        """
        return _spec.f90wrap_inputlist__get__wmanual()
    
    @wmanual.setter
    def wmanual(self, wmanual):
        _spec.f90wrap_inputlist__set__wmanual(wmanual)
    
    @property
    def wrzaxis(self):
        """
        Element wrzaxis ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 169
        
        """
        return _spec.f90wrap_inputlist__get__wrzaxis()
    
    @wrzaxis.setter
    def wrzaxis(self, wrzaxis):
        _spec.f90wrap_inputlist__set__wrzaxis(wrzaxis)
    
    @property
    def wpackxi(self):
        """
        Element wpackxi ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 170
        
        """
        return _spec.f90wrap_inputlist__get__wpackxi()
    
    @wpackxi.setter
    def wpackxi(self, wpackxi):
        _spec.f90wrap_inputlist__set__wpackxi(wpackxi)
    
    @property
    def wvolume(self):
        """
        Element wvolume ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 171
        
        """
        return _spec.f90wrap_inputlist__get__wvolume()
    
    @wvolume.setter
    def wvolume(self, wvolume):
        _spec.f90wrap_inputlist__set__wvolume(wvolume)
    
    @property
    def wcoords(self):
        """
        Element wcoords ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 172
        
        """
        return _spec.f90wrap_inputlist__get__wcoords()
    
    @wcoords.setter
    def wcoords(self, wcoords):
        _spec.f90wrap_inputlist__set__wcoords(wcoords)
    
    @property
    def wbasefn(self):
        """
        Element wbasefn ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 173
        
        """
        return _spec.f90wrap_inputlist__get__wbasefn()
    
    @wbasefn.setter
    def wbasefn(self, wbasefn):
        _spec.f90wrap_inputlist__set__wbasefn(wbasefn)
    
    @property
    def wmemory(self):
        """
        Element wmemory ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 174
        
        """
        return _spec.f90wrap_inputlist__get__wmemory()
    
    @wmemory.setter
    def wmemory(self, wmemory):
        _spec.f90wrap_inputlist__set__wmemory(wmemory)
    
    @property
    def wmetrix(self):
        """
        Element wmetrix ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 175
        
        """
        return _spec.f90wrap_inputlist__get__wmetrix()
    
    @wmetrix.setter
    def wmetrix(self, wmetrix):
        _spec.f90wrap_inputlist__set__wmetrix(wmetrix)
    
    @property
    def wma00aa(self):
        """
        Element wma00aa ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 176
        
        """
        return _spec.f90wrap_inputlist__get__wma00aa()
    
    @wma00aa.setter
    def wma00aa(self, wma00aa):
        _spec.f90wrap_inputlist__set__wma00aa(wma00aa)
    
    @property
    def wmatrix(self):
        """
        Element wmatrix ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 177
        
        """
        return _spec.f90wrap_inputlist__get__wmatrix()
    
    @wmatrix.setter
    def wmatrix(self, wmatrix):
        _spec.f90wrap_inputlist__set__wmatrix(wmatrix)
    
    @property
    def wspsmat(self):
        """
        Element wspsmat ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 178
        
        """
        return _spec.f90wrap_inputlist__get__wspsmat()
    
    @wspsmat.setter
    def wspsmat(self, wspsmat):
        _spec.f90wrap_inputlist__set__wspsmat(wspsmat)
    
    @property
    def wspsint(self):
        """
        Element wspsint ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 179
        
        """
        return _spec.f90wrap_inputlist__get__wspsint()
    
    @wspsint.setter
    def wspsint(self, wspsint):
        _spec.f90wrap_inputlist__set__wspsint(wspsint)
    
    @property
    def wmp00ac(self):
        """
        Element wmp00ac ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 180
        
        """
        return _spec.f90wrap_inputlist__get__wmp00ac()
    
    @wmp00ac.setter
    def wmp00ac(self, wmp00ac):
        _spec.f90wrap_inputlist__set__wmp00ac(wmp00ac)
    
    @property
    def wma02aa(self):
        """
        Element wma02aa ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 181
        
        """
        return _spec.f90wrap_inputlist__get__wma02aa()
    
    @wma02aa.setter
    def wma02aa(self, wma02aa):
        _spec.f90wrap_inputlist__set__wma02aa(wma02aa)
    
    @property
    def wpackab(self):
        """
        Element wpackab ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 182
        
        """
        return _spec.f90wrap_inputlist__get__wpackab()
    
    @wpackab.setter
    def wpackab(self, wpackab):
        _spec.f90wrap_inputlist__set__wpackab(wpackab)
    
    @property
    def wtr00ab(self):
        """
        Element wtr00ab ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 183
        
        """
        return _spec.f90wrap_inputlist__get__wtr00ab()
    
    @wtr00ab.setter
    def wtr00ab(self, wtr00ab):
        _spec.f90wrap_inputlist__set__wtr00ab(wtr00ab)
    
    @property
    def wcurent(self):
        """
        Element wcurent ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 184
        
        """
        return _spec.f90wrap_inputlist__get__wcurent()
    
    @wcurent.setter
    def wcurent(self, wcurent):
        _spec.f90wrap_inputlist__set__wcurent(wcurent)
    
    @property
    def wdf00ab(self):
        """
        Element wdf00ab ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 185
        
        """
        return _spec.f90wrap_inputlist__get__wdf00ab()
    
    @wdf00ab.setter
    def wdf00ab(self, wdf00ab):
        _spec.f90wrap_inputlist__set__wdf00ab(wdf00ab)
    
    @property
    def wlforce(self):
        """
        Element wlforce ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 186
        
        """
        return _spec.f90wrap_inputlist__get__wlforce()
    
    @wlforce.setter
    def wlforce(self, wlforce):
        _spec.f90wrap_inputlist__set__wlforce(wlforce)
    
    @property
    def wintghs(self):
        """
        Element wintghs ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 187
        
        """
        return _spec.f90wrap_inputlist__get__wintghs()
    
    @wintghs.setter
    def wintghs(self, wintghs):
        _spec.f90wrap_inputlist__set__wintghs(wintghs)
    
    @property
    def wmtrxhs(self):
        """
        Element wmtrxhs ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 188
        
        """
        return _spec.f90wrap_inputlist__get__wmtrxhs()
    
    @wmtrxhs.setter
    def wmtrxhs(self, wmtrxhs):
        _spec.f90wrap_inputlist__set__wmtrxhs(wmtrxhs)
    
    @property
    def wlbpol(self):
        """
        Element wlbpol ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 189
        
        """
        return _spec.f90wrap_inputlist__get__wlbpol()
    
    @wlbpol.setter
    def wlbpol(self, wlbpol):
        _spec.f90wrap_inputlist__set__wlbpol(wlbpol)
    
    @property
    def wbrcast(self):
        """
        Element wbrcast ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 190
        
        """
        return _spec.f90wrap_inputlist__get__wbrcast()
    
    @wbrcast.setter
    def wbrcast(self, wbrcast):
        _spec.f90wrap_inputlist__set__wbrcast(wbrcast)
    
    @property
    def wdfp100(self):
        """
        Element wdfp100 ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 191
        
        """
        return _spec.f90wrap_inputlist__get__wdfp100()
    
    @wdfp100.setter
    def wdfp100(self, wdfp100):
        _spec.f90wrap_inputlist__set__wdfp100(wdfp100)
    
    @property
    def wdfp200(self):
        """
        Element wdfp200 ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 192
        
        """
        return _spec.f90wrap_inputlist__get__wdfp200()
    
    @wdfp200.setter
    def wdfp200(self, wdfp200):
        _spec.f90wrap_inputlist__set__wdfp200(wdfp200)
    
    @property
    def wdforce(self):
        """
        Element wdforce ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 193
        
        """
        return _spec.f90wrap_inputlist__get__wdforce()
    
    @wdforce.setter
    def wdforce(self, wdforce):
        _spec.f90wrap_inputlist__set__wdforce(wdforce)
    
    @property
    def wnewton(self):
        """
        Element wnewton ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 194
        
        """
        return _spec.f90wrap_inputlist__get__wnewton()
    
    @wnewton.setter
    def wnewton(self, wnewton):
        _spec.f90wrap_inputlist__set__wnewton(wnewton)
    
    @property
    def wcasing(self):
        """
        Element wcasing ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 195
        
        """
        return _spec.f90wrap_inputlist__get__wcasing()
    
    @wcasing.setter
    def wcasing(self, wcasing):
        _spec.f90wrap_inputlist__set__wcasing(wcasing)
    
    @property
    def wbnorml(self):
        """
        Element wbnorml ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 196
        
        """
        return _spec.f90wrap_inputlist__get__wbnorml()
    
    @wbnorml.setter
    def wbnorml(self, wbnorml):
        _spec.f90wrap_inputlist__set__wbnorml(wbnorml)
    
    @property
    def wjo00aa(self):
        """
        Element wjo00aa ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 197
        
        """
        return _spec.f90wrap_inputlist__get__wjo00aa()
    
    @wjo00aa.setter
    def wjo00aa(self, wjo00aa):
        _spec.f90wrap_inputlist__set__wjo00aa(wjo00aa)
    
    @property
    def wpp00aa(self):
        """
        Element wpp00aa ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 198
        
        """
        return _spec.f90wrap_inputlist__get__wpp00aa()
    
    @wpp00aa.setter
    def wpp00aa(self, wpp00aa):
        _spec.f90wrap_inputlist__set__wpp00aa(wpp00aa)
    
    @property
    def wpp00ab(self):
        """
        Element wpp00ab ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 199
        
        """
        return _spec.f90wrap_inputlist__get__wpp00ab()
    
    @wpp00ab.setter
    def wpp00ab(self, wpp00ab):
        _spec.f90wrap_inputlist__set__wpp00ab(wpp00ab)
    
    @property
    def wbfield(self):
        """
        Element wbfield ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 200
        
        """
        return _spec.f90wrap_inputlist__get__wbfield()
    
    @wbfield.setter
    def wbfield(self, wbfield):
        _spec.f90wrap_inputlist__set__wbfield(wbfield)
    
    @property
    def wstzxyz(self):
        """
        Element wstzxyz ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 201
        
        """
        return _spec.f90wrap_inputlist__get__wstzxyz()
    
    @wstzxyz.setter
    def wstzxyz(self, wstzxyz):
        _spec.f90wrap_inputlist__set__wstzxyz(wstzxyz)
    
    @property
    def whesian(self):
        """
        Element whesian ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 202
        
        """
        return _spec.f90wrap_inputlist__get__whesian()
    
    @whesian.setter
    def whesian(self, whesian):
        _spec.f90wrap_inputlist__set__whesian(whesian)
    
    @property
    def wra00aa(self):
        """
        Element wra00aa ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 203
        
        """
        return _spec.f90wrap_inputlist__get__wra00aa()
    
    @wra00aa.setter
    def wra00aa(self, wra00aa):
        _spec.f90wrap_inputlist__set__wra00aa(wra00aa)
    
    @property
    def wnumrec(self):
        """
        Element wnumrec ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 204
        
        """
        return _spec.f90wrap_inputlist__get__wnumrec()
    
    @wnumrec.setter
    def wnumrec(self, wnumrec):
        _spec.f90wrap_inputlist__set__wnumrec(wnumrec)
    
    @property
    def wxspech(self):
        """
        Element wxspech ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 205
        
        """
        return _spec.f90wrap_inputlist__get__wxspech()
    
    @wxspech.setter
    def wxspech(self, wxspech):
        _spec.f90wrap_inputlist__set__wxspech(wxspech)
    
    @property
    def wbuild_vector_potential(self):
        """
        Element wbuild_vector_potential ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 207
        
        """
        return _spec.f90wrap_inputlist__get__wbuild_vector_potential()
    
    @wbuild_vector_potential.setter
    def wbuild_vector_potential(self, wbuild_vector_potential):
        _spec.f90wrap_inputlist__set__wbuild_vector_potential(wbuild_vector_potential)
    
    @property
    def wreadin(self):
        """
        Element wreadin ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 208
        
        """
        return _spec.f90wrap_inputlist__get__wreadin()
    
    @wreadin.setter
    def wreadin(self, wreadin):
        _spec.f90wrap_inputlist__set__wreadin(wreadin)
    
    @property
    def wwritin(self):
        """
        Element wwritin ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 209
        
        """
        return _spec.f90wrap_inputlist__get__wwritin()
    
    @wwritin.setter
    def wwritin(self, wwritin):
        _spec.f90wrap_inputlist__set__wwritin(wwritin)
    
    @property
    def wwrtend(self):
        """
        Element wwrtend ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 210
        
        """
        return _spec.f90wrap_inputlist__get__wwrtend()
    
    @wwrtend.setter
    def wwrtend(self, wwrtend):
        _spec.f90wrap_inputlist__set__wwrtend(wwrtend)
    
    @property
    def wmacros(self):
        """
        Element wmacros ftype=logical pytype=bool
        
        
        Defined at inputlist.fpp line 211
        
        """
        return _spec.f90wrap_inputlist__get__wmacros()
    
    @wmacros.setter
    def wmacros(self, wmacros):
        _spec.f90wrap_inputlist__set__wmacros(wmacros)
    
    def __str__(self):
        ret = ['<inputlist>{\n']
        ret.append('    mnvol : ')
        ret.append(repr(self.mnvol))
        ret.append(',\n    mmpol : ')
        ret.append(repr(self.mmpol))
        ret.append(',\n    mntor : ')
        ret.append(repr(self.mntor))
        ret.append(',\n    igeometry : ')
        ret.append(repr(self.igeometry))
        ret.append(',\n    istellsym : ')
        ret.append(repr(self.istellsym))
        ret.append(',\n    lfreebound : ')
        ret.append(repr(self.lfreebound))
        ret.append(',\n    phiedge : ')
        ret.append(repr(self.phiedge))
        ret.append(',\n    curtor : ')
        ret.append(repr(self.curtor))
        ret.append(',\n    curpol : ')
        ret.append(repr(self.curpol))
        ret.append(',\n    gamma : ')
        ret.append(repr(self.gamma))
        ret.append(',\n    nfp : ')
        ret.append(repr(self.nfp))
        ret.append(',\n    nvol : ')
        ret.append(repr(self.nvol))
        ret.append(',\n    mpol : ')
        ret.append(repr(self.mpol))
        ret.append(',\n    ntor : ')
        ret.append(repr(self.ntor))
        ret.append(',\n    lrad : ')
        ret.append(repr(self.lrad))
        ret.append(',\n    lconstraint : ')
        ret.append(repr(self.lconstraint))
        ret.append(',\n    tflux : ')
        ret.append(repr(self.tflux))
        ret.append(',\n    pflux : ')
        ret.append(repr(self.pflux))
        ret.append(',\n    helicity : ')
        ret.append(repr(self.helicity))
        ret.append(',\n    pscale : ')
        ret.append(repr(self.pscale))
        ret.append(',\n    pressure : ')
        ret.append(repr(self.pressure))
        ret.append(',\n    ladiabatic : ')
        ret.append(repr(self.ladiabatic))
        ret.append(',\n    adiabatic : ')
        ret.append(repr(self.adiabatic))
        ret.append(',\n    mu : ')
        ret.append(repr(self.mu))
        ret.append(',\n    ivolume : ')
        ret.append(repr(self.ivolume))
        ret.append(',\n    isurf : ')
        ret.append(repr(self.isurf))
        ret.append(',\n    pl : ')
        ret.append(repr(self.pl))
        ret.append(',\n    ql : ')
        ret.append(repr(self.ql))
        ret.append(',\n    pr : ')
        ret.append(repr(self.pr))
        ret.append(',\n    qr : ')
        ret.append(repr(self.qr))
        ret.append(',\n    iota : ')
        ret.append(repr(self.iota))
        ret.append(',\n    lp : ')
        ret.append(repr(self.lp))
        ret.append(',\n    lq : ')
        ret.append(repr(self.lq))
        ret.append(',\n    rp : ')
        ret.append(repr(self.rp))
        ret.append(',\n    rq : ')
        ret.append(repr(self.rq))
        ret.append(',\n    oita : ')
        ret.append(repr(self.oita))
        ret.append(',\n    rpol : ')
        ret.append(repr(self.rpol))
        ret.append(',\n    rtor : ')
        ret.append(repr(self.rtor))
        ret.append(',\n    rac : ')
        ret.append(repr(self.rac))
        ret.append(',\n    zas : ')
        ret.append(repr(self.zas))
        ret.append(',\n    ras : ')
        ret.append(repr(self.ras))
        ret.append(',\n    zac : ')
        ret.append(repr(self.zac))
        ret.append(',\n    rbc : ')
        ret.append(repr(self.rbc))
        ret.append(',\n    zbs : ')
        ret.append(repr(self.zbs))
        ret.append(',\n    rbs : ')
        ret.append(repr(self.rbs))
        ret.append(',\n    zbc : ')
        ret.append(repr(self.zbc))
        ret.append(',\n    rwc : ')
        ret.append(repr(self.rwc))
        ret.append(',\n    zws : ')
        ret.append(repr(self.zws))
        ret.append(',\n    rws : ')
        ret.append(repr(self.rws))
        ret.append(',\n    zwc : ')
        ret.append(repr(self.zwc))
        ret.append(',\n    vns : ')
        ret.append(repr(self.vns))
        ret.append(',\n    bns : ')
        ret.append(repr(self.bns))
        ret.append(',\n    vnc : ')
        ret.append(repr(self.vnc))
        ret.append(',\n    bnc : ')
        ret.append(repr(self.bnc))
        ret.append(',\n    mupftol : ')
        ret.append(repr(self.mupftol))
        ret.append(',\n    mupfits : ')
        ret.append(repr(self.mupfits))
        ret.append(',\n    lreflect : ')
        ret.append(repr(self.lreflect))
        ret.append(',\n    linitialize : ')
        ret.append(repr(self.linitialize))
        ret.append(',\n    lautoinitbn : ')
        ret.append(repr(self.lautoinitbn))
        ret.append(',\n    lzerovac : ')
        ret.append(repr(self.lzerovac))
        ret.append(',\n    ndiscrete : ')
        ret.append(repr(self.ndiscrete))
        ret.append(',\n    nquad : ')
        ret.append(repr(self.nquad))
        ret.append(',\n    impol : ')
        ret.append(repr(self.impol))
        ret.append(',\n    intor : ')
        ret.append(repr(self.intor))
        ret.append(',\n    lsparse : ')
        ret.append(repr(self.lsparse))
        ret.append(',\n    lsvdiota : ')
        ret.append(repr(self.lsvdiota))
        ret.append(',\n    imethod : ')
        ret.append(repr(self.imethod))
        ret.append(',\n    iorder : ')
        ret.append(repr(self.iorder))
        ret.append(',\n    iprecon : ')
        ret.append(repr(self.iprecon))
        ret.append(',\n    iotatol : ')
        ret.append(repr(self.iotatol))
        ret.append(',\n    lextrap : ')
        ret.append(repr(self.lextrap))
        ret.append(',\n    mregular : ')
        ret.append(repr(self.mregular))
        ret.append(',\n    lrzaxis : ')
        ret.append(repr(self.lrzaxis))
        ret.append(',\n    ntoraxis : ')
        ret.append(repr(self.ntoraxis))
        ret.append(',\n    lbeltrami : ')
        ret.append(repr(self.lbeltrami))
        ret.append(',\n    linitgues : ')
        ret.append(repr(self.linitgues))
        ret.append(',\n    lposdef : ')
        ret.append(repr(self.lposdef))
        ret.append(',\n    maxrndgues : ')
        ret.append(repr(self.maxrndgues))
        ret.append(',\n    lmatsolver : ')
        ret.append(repr(self.lmatsolver))
        ret.append(',\n    nitergmres : ')
        ret.append(repr(self.nitergmres))
        ret.append(',\n    epsgmres : ')
        ret.append(repr(self.epsgmres))
        ret.append(',\n    lgmresprec : ')
        ret.append(repr(self.lgmresprec))
        ret.append(',\n    epsilu : ')
        ret.append(repr(self.epsilu))
        ret.append(',\n    lfindzero : ')
        ret.append(repr(self.lfindzero))
        ret.append(',\n    escale : ')
        ret.append(repr(self.escale))
        ret.append(',\n    opsilon : ')
        ret.append(repr(self.opsilon))
        ret.append(',\n    pcondense : ')
        ret.append(repr(self.pcondense))
        ret.append(',\n    epsilon : ')
        ret.append(repr(self.epsilon))
        ret.append(',\n    wpoloidal : ')
        ret.append(repr(self.wpoloidal))
        ret.append(',\n    upsilon : ')
        ret.append(repr(self.upsilon))
        ret.append(',\n    forcetol : ')
        ret.append(repr(self.forcetol))
        ret.append(',\n    c05xmax : ')
        ret.append(repr(self.c05xmax))
        ret.append(',\n    c05xtol : ')
        ret.append(repr(self.c05xtol))
        ret.append(',\n    c05factor : ')
        ret.append(repr(self.c05factor))
        ret.append(',\n    lreadgf : ')
        ret.append(repr(self.lreadgf))
        ret.append(',\n    mfreeits : ')
        ret.append(repr(self.mfreeits))
        ret.append(',\n    bnstol : ')
        ret.append(repr(self.bnstol))
        ret.append(',\n    bnsblend : ')
        ret.append(repr(self.bnsblend))
        ret.append(',\n    gbntol : ')
        ret.append(repr(self.gbntol))
        ret.append(',\n    gbnbld : ')
        ret.append(repr(self.gbnbld))
        ret.append(',\n    vcasingeps : ')
        ret.append(repr(self.vcasingeps))
        ret.append(',\n    vcasingtol : ')
        ret.append(repr(self.vcasingtol))
        ret.append(',\n    vcasingits : ')
        ret.append(repr(self.vcasingits))
        ret.append(',\n    vcasingper : ')
        ret.append(repr(self.vcasingper))
        ret.append(',\n    mcasingcal : ')
        ret.append(repr(self.mcasingcal))
        ret.append(',\n    odetol : ')
        ret.append(repr(self.odetol))
        ret.append(',\n    absreq : ')
        ret.append(repr(self.absreq))
        ret.append(',\n    relreq : ')
        ret.append(repr(self.relreq))
        ret.append(',\n    absacc : ')
        ret.append(repr(self.absacc))
        ret.append(',\n    epsr : ')
        ret.append(repr(self.epsr))
        ret.append(',\n    nppts : ')
        ret.append(repr(self.nppts))
        ret.append(',\n    ppts : ')
        ret.append(repr(self.ppts))
        ret.append(',\n    nptrj : ')
        ret.append(repr(self.nptrj))
        ret.append(',\n    lhevalues : ')
        ret.append(repr(self.lhevalues))
        ret.append(',\n    lhevectors : ')
        ret.append(repr(self.lhevectors))
        ret.append(',\n    lhmatrix : ')
        ret.append(repr(self.lhmatrix))
        ret.append(',\n    lperturbed : ')
        ret.append(repr(self.lperturbed))
        ret.append(',\n    dpp : ')
        ret.append(repr(self.dpp))
        ret.append(',\n    dqq : ')
        ret.append(repr(self.dqq))
        ret.append(',\n    lerrortype : ')
        ret.append(repr(self.lerrortype))
        ret.append(',\n    ngrid : ')
        ret.append(repr(self.ngrid))
        ret.append(',\n    drz : ')
        ret.append(repr(self.drz))
        ret.append(',\n    lcheck : ')
        ret.append(repr(self.lcheck))
        ret.append(',\n    ltiming : ')
        ret.append(repr(self.ltiming))
        ret.append(',\n    fudge : ')
        ret.append(repr(self.fudge))
        ret.append(',\n    scaling : ')
        ret.append(repr(self.scaling))
        ret.append(',\n    wdcuhre : ')
        ret.append(repr(self.wdcuhre))
        ret.append(',\n    wminpack : ')
        ret.append(repr(self.wminpack))
        ret.append(',\n    wiqpack : ')
        ret.append(repr(self.wiqpack))
        ret.append(',\n    wrksuite : ')
        ret.append(repr(self.wrksuite))
        ret.append(',\n    wi1mach : ')
        ret.append(repr(self.wi1mach))
        ret.append(',\n    wd1mach : ')
        ret.append(repr(self.wd1mach))
        ret.append(',\n    wilut : ')
        ret.append(repr(self.wilut))
        ret.append(',\n    witers : ')
        ret.append(repr(self.witers))
        ret.append(',\n    winputlist : ')
        ret.append(repr(self.winputlist))
        ret.append(',\n    wglobal : ')
        ret.append(repr(self.wglobal))
        ret.append(',\n    wsphdf5 : ')
        ret.append(repr(self.wsphdf5))
        ret.append(',\n    wpreset : ')
        ret.append(repr(self.wpreset))
        ret.append(',\n    wmanual : ')
        ret.append(repr(self.wmanual))
        ret.append(',\n    wrzaxis : ')
        ret.append(repr(self.wrzaxis))
        ret.append(',\n    wpackxi : ')
        ret.append(repr(self.wpackxi))
        ret.append(',\n    wvolume : ')
        ret.append(repr(self.wvolume))
        ret.append(',\n    wcoords : ')
        ret.append(repr(self.wcoords))
        ret.append(',\n    wbasefn : ')
        ret.append(repr(self.wbasefn))
        ret.append(',\n    wmemory : ')
        ret.append(repr(self.wmemory))
        ret.append(',\n    wmetrix : ')
        ret.append(repr(self.wmetrix))
        ret.append(',\n    wma00aa : ')
        ret.append(repr(self.wma00aa))
        ret.append(',\n    wmatrix : ')
        ret.append(repr(self.wmatrix))
        ret.append(',\n    wspsmat : ')
        ret.append(repr(self.wspsmat))
        ret.append(',\n    wspsint : ')
        ret.append(repr(self.wspsint))
        ret.append(',\n    wmp00ac : ')
        ret.append(repr(self.wmp00ac))
        ret.append(',\n    wma02aa : ')
        ret.append(repr(self.wma02aa))
        ret.append(',\n    wpackab : ')
        ret.append(repr(self.wpackab))
        ret.append(',\n    wtr00ab : ')
        ret.append(repr(self.wtr00ab))
        ret.append(',\n    wcurent : ')
        ret.append(repr(self.wcurent))
        ret.append(',\n    wdf00ab : ')
        ret.append(repr(self.wdf00ab))
        ret.append(',\n    wlforce : ')
        ret.append(repr(self.wlforce))
        ret.append(',\n    wintghs : ')
        ret.append(repr(self.wintghs))
        ret.append(',\n    wmtrxhs : ')
        ret.append(repr(self.wmtrxhs))
        ret.append(',\n    wlbpol : ')
        ret.append(repr(self.wlbpol))
        ret.append(',\n    wbrcast : ')
        ret.append(repr(self.wbrcast))
        ret.append(',\n    wdfp100 : ')
        ret.append(repr(self.wdfp100))
        ret.append(',\n    wdfp200 : ')
        ret.append(repr(self.wdfp200))
        ret.append(',\n    wdforce : ')
        ret.append(repr(self.wdforce))
        ret.append(',\n    wnewton : ')
        ret.append(repr(self.wnewton))
        ret.append(',\n    wcasing : ')
        ret.append(repr(self.wcasing))
        ret.append(',\n    wbnorml : ')
        ret.append(repr(self.wbnorml))
        ret.append(',\n    wjo00aa : ')
        ret.append(repr(self.wjo00aa))
        ret.append(',\n    wpp00aa : ')
        ret.append(repr(self.wpp00aa))
        ret.append(',\n    wpp00ab : ')
        ret.append(repr(self.wpp00ab))
        ret.append(',\n    wbfield : ')
        ret.append(repr(self.wbfield))
        ret.append(',\n    wstzxyz : ')
        ret.append(repr(self.wstzxyz))
        ret.append(',\n    whesian : ')
        ret.append(repr(self.whesian))
        ret.append(',\n    wra00aa : ')
        ret.append(repr(self.wra00aa))
        ret.append(',\n    wnumrec : ')
        ret.append(repr(self.wnumrec))
        ret.append(',\n    wxspech : ')
        ret.append(repr(self.wxspech))
        ret.append(',\n    wbuild_vector_potential : ')
        ret.append(repr(self.wbuild_vector_potential))
        ret.append(',\n    wreadin : ')
        ret.append(repr(self.wreadin))
        ret.append(',\n    wwritin : ')
        ret.append(repr(self.wwritin))
        ret.append(',\n    wwrtend : ')
        ret.append(repr(self.wwrtend))
        ret.append(',\n    wmacros : ')
        ret.append(repr(self.wmacros))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

inputlist = Inputlist()

class Constants(f90wrap.runtime.FortranModule):
    """
    Module constants
    
    
    Defined at global.fpp lines 24-51
    
    """
    @property
    def zero(self):
        """
        Element zero ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 26
        
        """
        return _spec.f90wrap_constants__get__zero()
    
    @property
    def one(self):
        """
        Element one ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 27
        
        """
        return _spec.f90wrap_constants__get__one()
    
    @property
    def two(self):
        """
        Element two ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 28
        
        """
        return _spec.f90wrap_constants__get__two()
    
    @property
    def three(self):
        """
        Element three ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 29
        
        """
        return _spec.f90wrap_constants__get__three()
    
    @property
    def four(self):
        """
        Element four ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 30
        
        """
        return _spec.f90wrap_constants__get__four()
    
    @property
    def five(self):
        """
        Element five ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 31
        
        """
        return _spec.f90wrap_constants__get__five()
    
    @property
    def six(self):
        """
        Element six ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 32
        
        """
        return _spec.f90wrap_constants__get__six()
    
    @property
    def seven(self):
        """
        Element seven ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 33
        
        """
        return _spec.f90wrap_constants__get__seven()
    
    @property
    def eight(self):
        """
        Element eight ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 34
        
        """
        return _spec.f90wrap_constants__get__eight()
    
    @property
    def nine(self):
        """
        Element nine ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 35
        
        """
        return _spec.f90wrap_constants__get__nine()
    
    @property
    def ten(self):
        """
        Element ten ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 36
        
        """
        return _spec.f90wrap_constants__get__ten()
    
    @property
    def eleven(self):
        """
        Element eleven ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 37
        
        """
        return _spec.f90wrap_constants__get__eleven()
    
    @property
    def twelve(self):
        """
        Element twelve ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 38
        
        """
        return _spec.f90wrap_constants__get__twelve()
    
    @property
    def hundred(self):
        """
        Element hundred ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 39
        
        """
        return _spec.f90wrap_constants__get__hundred()
    
    @property
    def thousand(self):
        """
        Element thousand ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 40
        
        """
        return _spec.f90wrap_constants__get__thousand()
    
    @property
    def half(self):
        """
        Element half ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 41
        
        """
        return _spec.f90wrap_constants__get__half()
    
    @property
    def third(self):
        """
        Element third ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 42
        
        """
        return _spec.f90wrap_constants__get__third()
    
    @property
    def quart(self):
        """
        Element quart ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 43
        
        """
        return _spec.f90wrap_constants__get__quart()
    
    @property
    def fifth(self):
        """
        Element fifth ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 44
        
        """
        return _spec.f90wrap_constants__get__fifth()
    
    @property
    def sixth(self):
        """
        Element sixth ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 45
        
        """
        return _spec.f90wrap_constants__get__sixth()
    
    @property
    def pi2(self):
        """
        Element pi2 ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 46
        
        """
        return _spec.f90wrap_constants__get__pi2()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 47
        
        """
        return _spec.f90wrap_constants__get__pi()
    
    @property
    def mu0(self):
        """
        Element mu0 ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 48
        
        """
        return _spec.f90wrap_constants__get__mu0()
    
    @property
    def goldenmean(self):
        """
        Element goldenmean ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 49
        
        """
        return _spec.f90wrap_constants__get__goldenmean()
    
    @property
    def version(self):
        """
        Element version ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 51
        
        """
        return _spec.f90wrap_constants__get__version()
    
    def __str__(self):
        ret = ['<constants>{\n']
        ret.append('    zero : ')
        ret.append(repr(self.zero))
        ret.append(',\n    one : ')
        ret.append(repr(self.one))
        ret.append(',\n    two : ')
        ret.append(repr(self.two))
        ret.append(',\n    three : ')
        ret.append(repr(self.three))
        ret.append(',\n    four : ')
        ret.append(repr(self.four))
        ret.append(',\n    five : ')
        ret.append(repr(self.five))
        ret.append(',\n    six : ')
        ret.append(repr(self.six))
        ret.append(',\n    seven : ')
        ret.append(repr(self.seven))
        ret.append(',\n    eight : ')
        ret.append(repr(self.eight))
        ret.append(',\n    nine : ')
        ret.append(repr(self.nine))
        ret.append(',\n    ten : ')
        ret.append(repr(self.ten))
        ret.append(',\n    eleven : ')
        ret.append(repr(self.eleven))
        ret.append(',\n    twelve : ')
        ret.append(repr(self.twelve))
        ret.append(',\n    hundred : ')
        ret.append(repr(self.hundred))
        ret.append(',\n    thousand : ')
        ret.append(repr(self.thousand))
        ret.append(',\n    half : ')
        ret.append(repr(self.half))
        ret.append(',\n    third : ')
        ret.append(repr(self.third))
        ret.append(',\n    quart : ')
        ret.append(repr(self.quart))
        ret.append(',\n    fifth : ')
        ret.append(repr(self.fifth))
        ret.append(',\n    sixth : ')
        ret.append(repr(self.sixth))
        ret.append(',\n    pi2 : ')
        ret.append(repr(self.pi2))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    mu0 : ')
        ret.append(repr(self.mu0))
        ret.append(',\n    goldenmean : ')
        ret.append(repr(self.goldenmean))
        ret.append(',\n    version : ')
        ret.append(repr(self.version))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

constants = Constants()

class Numerical(f90wrap.runtime.FortranModule):
    """
    Module numerical
    
    
    Defined at global.fpp lines 55-65
    
    """
    @property
    def machprec(self):
        """
        Element machprec ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 61
        
        """
        return _spec.f90wrap_numerical__get__machprec()
    
    @property
    def vsmall(self):
        """
        Element vsmall ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 62
        
        """
        return _spec.f90wrap_numerical__get__vsmall()
    
    @property
    def small(self):
        """
        Element small ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 63
        
        """
        return _spec.f90wrap_numerical__get__small()
    
    @property
    def sqrtmachprec(self):
        """
        Element sqrtmachprec ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 64
        
        """
        return _spec.f90wrap_numerical__get__sqrtmachprec()
    
    @property
    def logtolerance(self):
        """
        Element logtolerance ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 65
        
        """
        return _spec.f90wrap_numerical__get__logtolerance()
    
    def __str__(self):
        ret = ['<numerical>{\n']
        ret.append('    machprec : ')
        ret.append(repr(self.machprec))
        ret.append(',\n    vsmall : ')
        ret.append(repr(self.vsmall))
        ret.append(',\n    small : ')
        ret.append(repr(self.small))
        ret.append(',\n    sqrtmachprec : ')
        ret.append(repr(self.sqrtmachprec))
        ret.append(',\n    logtolerance : ')
        ret.append(repr(self.logtolerance))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

numerical = Numerical()

class Fileunits(f90wrap.runtime.FortranModule):
    """
    Module fileunits
    
    
    Defined at global.fpp lines 68-97
    
    """
    @staticmethod
    def mute(action):
        """
        mute(action)
        
        
        Defined at global.fpp lines 81-96
        
        Parameters
        ----------
        action : int
        
        """
        _spec.f90wrap_mute(action=action)
    
    @property
    def iunit(self):
        """
        Element iunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 70
        
        """
        return _spec.f90wrap_fileunits__get__iunit()
    
    @iunit.setter
    def iunit(self, iunit):
        _spec.f90wrap_fileunits__set__iunit(iunit)
    
    @property
    def ounit(self):
        """
        Element ounit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 71
        
        """
        return _spec.f90wrap_fileunits__get__ounit()
    
    @ounit.setter
    def ounit(self, ounit):
        _spec.f90wrap_fileunits__set__ounit(ounit)
    
    @property
    def gunit(self):
        """
        Element gunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 72
        
        """
        return _spec.f90wrap_fileunits__get__gunit()
    
    @gunit.setter
    def gunit(self, gunit):
        _spec.f90wrap_fileunits__set__gunit(gunit)
    
    @property
    def aunit(self):
        """
        Element aunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 73
        
        """
        return _spec.f90wrap_fileunits__get__aunit()
    
    @aunit.setter
    def aunit(self, aunit):
        _spec.f90wrap_fileunits__set__aunit(aunit)
    
    @property
    def dunit(self):
        """
        Element dunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 74
        
        """
        return _spec.f90wrap_fileunits__get__dunit()
    
    @dunit.setter
    def dunit(self, dunit):
        _spec.f90wrap_fileunits__set__dunit(dunit)
    
    @property
    def hunit(self):
        """
        Element hunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 75
        
        """
        return _spec.f90wrap_fileunits__get__hunit()
    
    @hunit.setter
    def hunit(self, hunit):
        _spec.f90wrap_fileunits__set__hunit(hunit)
    
    @property
    def munit(self):
        """
        Element munit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 76
        
        """
        return _spec.f90wrap_fileunits__get__munit()
    
    @munit.setter
    def munit(self, munit):
        _spec.f90wrap_fileunits__set__munit(munit)
    
    @property
    def lunit(self):
        """
        Element lunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 77
        
        """
        return _spec.f90wrap_fileunits__get__lunit()
    
    @lunit.setter
    def lunit(self, lunit):
        _spec.f90wrap_fileunits__set__lunit(lunit)
    
    @property
    def vunit(self):
        """
        Element vunit ftype=integer  pytype=int
        
        
        Defined at global.fpp line 78
        
        """
        return _spec.f90wrap_fileunits__get__vunit()
    
    @vunit.setter
    def vunit(self, vunit):
        _spec.f90wrap_fileunits__set__vunit(vunit)
    
    def __str__(self):
        ret = ['<fileunits>{\n']
        ret.append('    iunit : ')
        ret.append(repr(self.iunit))
        ret.append(',\n    ounit : ')
        ret.append(repr(self.ounit))
        ret.append(',\n    gunit : ')
        ret.append(repr(self.gunit))
        ret.append(',\n    aunit : ')
        ret.append(repr(self.aunit))
        ret.append(',\n    dunit : ')
        ret.append(repr(self.dunit))
        ret.append(',\n    hunit : ')
        ret.append(repr(self.hunit))
        ret.append(',\n    munit : ')
        ret.append(repr(self.munit))
        ret.append(',\n    lunit : ')
        ret.append(repr(self.lunit))
        ret.append(',\n    vunit : ')
        ret.append(repr(self.vunit))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

fileunits = Fileunits()

class Cputiming(f90wrap.runtime.FortranModule):
    """
    Module cputiming
    
    
    Defined at global.fpp lines 100-154
    
    """
    @property
    def tdcuhre(self):
        """
        Element tdcuhre ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 101
        
        """
        return _spec.f90wrap_cputiming__get__tdcuhre()
    
    @tdcuhre.setter
    def tdcuhre(self, tdcuhre):
        _spec.f90wrap_cputiming__set__tdcuhre(tdcuhre)
    
    @property
    def dcuhret(self):
        """
        Element dcuhret ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 101
        
        """
        return _spec.f90wrap_cputiming__get__dcuhret()
    
    @dcuhret.setter
    def dcuhret(self, dcuhret):
        _spec.f90wrap_cputiming__set__dcuhret(dcuhret)
    
    @property
    def tminpack(self):
        """
        Element tminpack ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 102
        
        """
        return _spec.f90wrap_cputiming__get__tminpack()
    
    @tminpack.setter
    def tminpack(self, tminpack):
        _spec.f90wrap_cputiming__set__tminpack(tminpack)
    
    @property
    def minpackt(self):
        """
        Element minpackt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 102
        
        """
        return _spec.f90wrap_cputiming__get__minpackt()
    
    @minpackt.setter
    def minpackt(self, minpackt):
        _spec.f90wrap_cputiming__set__minpackt(minpackt)
    
    @property
    def tiqpack(self):
        """
        Element tiqpack ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 103
        
        """
        return _spec.f90wrap_cputiming__get__tiqpack()
    
    @tiqpack.setter
    def tiqpack(self, tiqpack):
        _spec.f90wrap_cputiming__set__tiqpack(tiqpack)
    
    @property
    def iqpackt(self):
        """
        Element iqpackt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 103
        
        """
        return _spec.f90wrap_cputiming__get__iqpackt()
    
    @iqpackt.setter
    def iqpackt(self, iqpackt):
        _spec.f90wrap_cputiming__set__iqpackt(iqpackt)
    
    @property
    def trksuite(self):
        """
        Element trksuite ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 104
        
        """
        return _spec.f90wrap_cputiming__get__trksuite()
    
    @trksuite.setter
    def trksuite(self, trksuite):
        _spec.f90wrap_cputiming__set__trksuite(trksuite)
    
    @property
    def rksuitet(self):
        """
        Element rksuitet ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 104
        
        """
        return _spec.f90wrap_cputiming__get__rksuitet()
    
    @rksuitet.setter
    def rksuitet(self, rksuitet):
        _spec.f90wrap_cputiming__set__rksuitet(rksuitet)
    
    @property
    def ti1mach(self):
        """
        Element ti1mach ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 105
        
        """
        return _spec.f90wrap_cputiming__get__ti1mach()
    
    @ti1mach.setter
    def ti1mach(self, ti1mach):
        _spec.f90wrap_cputiming__set__ti1mach(ti1mach)
    
    @property
    def i1macht(self):
        """
        Element i1macht ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 105
        
        """
        return _spec.f90wrap_cputiming__get__i1macht()
    
    @i1macht.setter
    def i1macht(self, i1macht):
        _spec.f90wrap_cputiming__set__i1macht(i1macht)
    
    @property
    def td1mach(self):
        """
        Element td1mach ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 106
        
        """
        return _spec.f90wrap_cputiming__get__td1mach()
    
    @td1mach.setter
    def td1mach(self, td1mach):
        _spec.f90wrap_cputiming__set__td1mach(td1mach)
    
    @property
    def d1macht(self):
        """
        Element d1macht ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 106
        
        """
        return _spec.f90wrap_cputiming__get__d1macht()
    
    @d1macht.setter
    def d1macht(self, d1macht):
        _spec.f90wrap_cputiming__set__d1macht(d1macht)
    
    @property
    def tilut(self):
        """
        Element tilut ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 107
        
        """
        return _spec.f90wrap_cputiming__get__tilut()
    
    @tilut.setter
    def tilut(self, tilut):
        _spec.f90wrap_cputiming__set__tilut(tilut)
    
    @property
    def ilutt(self):
        """
        Element ilutt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 107
        
        """
        return _spec.f90wrap_cputiming__get__ilutt()
    
    @ilutt.setter
    def ilutt(self, ilutt):
        _spec.f90wrap_cputiming__set__ilutt(ilutt)
    
    @property
    def titers(self):
        """
        Element titers ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 108
        
        """
        return _spec.f90wrap_cputiming__get__titers()
    
    @titers.setter
    def titers(self, titers):
        _spec.f90wrap_cputiming__set__titers(titers)
    
    @property
    def iterst(self):
        """
        Element iterst ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 108
        
        """
        return _spec.f90wrap_cputiming__get__iterst()
    
    @iterst.setter
    def iterst(self, iterst):
        _spec.f90wrap_cputiming__set__iterst(iterst)
    
    @property
    def tinputlist(self):
        """
        Element tinputlist ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 109
        
        """
        return _spec.f90wrap_cputiming__get__tinputlist()
    
    @tinputlist.setter
    def tinputlist(self, tinputlist):
        _spec.f90wrap_cputiming__set__tinputlist(tinputlist)
    
    @property
    def inputlistt(self):
        """
        Element inputlistt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 109
        
        """
        return _spec.f90wrap_cputiming__get__inputlistt()
    
    @inputlistt.setter
    def inputlistt(self, inputlistt):
        _spec.f90wrap_cputiming__set__inputlistt(inputlistt)
    
    @property
    def tglobal(self):
        """
        Element tglobal ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 110
        
        """
        return _spec.f90wrap_cputiming__get__tglobal()
    
    @tglobal.setter
    def tglobal(self, tglobal):
        _spec.f90wrap_cputiming__set__tglobal(tglobal)
    
    @property
    def globalt(self):
        """
        Element globalt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 110
        
        """
        return _spec.f90wrap_cputiming__get__globalt()
    
    @globalt.setter
    def globalt(self, globalt):
        _spec.f90wrap_cputiming__set__globalt(globalt)
    
    @property
    def tsphdf5(self):
        """
        Element tsphdf5 ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 111
        
        """
        return _spec.f90wrap_cputiming__get__tsphdf5()
    
    @tsphdf5.setter
    def tsphdf5(self, tsphdf5):
        _spec.f90wrap_cputiming__set__tsphdf5(tsphdf5)
    
    @property
    def sphdf5t(self):
        """
        Element sphdf5t ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 111
        
        """
        return _spec.f90wrap_cputiming__get__sphdf5t()
    
    @sphdf5t.setter
    def sphdf5t(self, sphdf5t):
        _spec.f90wrap_cputiming__set__sphdf5t(sphdf5t)
    
    @property
    def tpreset(self):
        """
        Element tpreset ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 112
        
        """
        return _spec.f90wrap_cputiming__get__tpreset()
    
    @tpreset.setter
    def tpreset(self, tpreset):
        _spec.f90wrap_cputiming__set__tpreset(tpreset)
    
    @property
    def presett(self):
        """
        Element presett ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 112
        
        """
        return _spec.f90wrap_cputiming__get__presett()
    
    @presett.setter
    def presett(self, presett):
        _spec.f90wrap_cputiming__set__presett(presett)
    
    @property
    def tmanual(self):
        """
        Element tmanual ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 113
        
        """
        return _spec.f90wrap_cputiming__get__tmanual()
    
    @tmanual.setter
    def tmanual(self, tmanual):
        _spec.f90wrap_cputiming__set__tmanual(tmanual)
    
    @property
    def manualt(self):
        """
        Element manualt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 113
        
        """
        return _spec.f90wrap_cputiming__get__manualt()
    
    @manualt.setter
    def manualt(self, manualt):
        _spec.f90wrap_cputiming__set__manualt(manualt)
    
    @property
    def trzaxis(self):
        """
        Element trzaxis ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 114
        
        """
        return _spec.f90wrap_cputiming__get__trzaxis()
    
    @trzaxis.setter
    def trzaxis(self, trzaxis):
        _spec.f90wrap_cputiming__set__trzaxis(trzaxis)
    
    @property
    def rzaxist(self):
        """
        Element rzaxist ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 114
        
        """
        return _spec.f90wrap_cputiming__get__rzaxist()
    
    @rzaxist.setter
    def rzaxist(self, rzaxist):
        _spec.f90wrap_cputiming__set__rzaxist(rzaxist)
    
    @property
    def tpackxi(self):
        """
        Element tpackxi ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 115
        
        """
        return _spec.f90wrap_cputiming__get__tpackxi()
    
    @tpackxi.setter
    def tpackxi(self, tpackxi):
        _spec.f90wrap_cputiming__set__tpackxi(tpackxi)
    
    @property
    def packxit(self):
        """
        Element packxit ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 115
        
        """
        return _spec.f90wrap_cputiming__get__packxit()
    
    @packxit.setter
    def packxit(self, packxit):
        _spec.f90wrap_cputiming__set__packxit(packxit)
    
    @property
    def tvolume(self):
        """
        Element tvolume ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 116
        
        """
        return _spec.f90wrap_cputiming__get__tvolume()
    
    @tvolume.setter
    def tvolume(self, tvolume):
        _spec.f90wrap_cputiming__set__tvolume(tvolume)
    
    @property
    def volumet(self):
        """
        Element volumet ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 116
        
        """
        return _spec.f90wrap_cputiming__get__volumet()
    
    @volumet.setter
    def volumet(self, volumet):
        _spec.f90wrap_cputiming__set__volumet(volumet)
    
    @property
    def tcoords(self):
        """
        Element tcoords ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 117
        
        """
        return _spec.f90wrap_cputiming__get__tcoords()
    
    @tcoords.setter
    def tcoords(self, tcoords):
        _spec.f90wrap_cputiming__set__tcoords(tcoords)
    
    @property
    def coordst(self):
        """
        Element coordst ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 117
        
        """
        return _spec.f90wrap_cputiming__get__coordst()
    
    @coordst.setter
    def coordst(self, coordst):
        _spec.f90wrap_cputiming__set__coordst(coordst)
    
    @property
    def tbasefn(self):
        """
        Element tbasefn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 118
        
        """
        return _spec.f90wrap_cputiming__get__tbasefn()
    
    @tbasefn.setter
    def tbasefn(self, tbasefn):
        _spec.f90wrap_cputiming__set__tbasefn(tbasefn)
    
    @property
    def basefnt(self):
        """
        Element basefnt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 118
        
        """
        return _spec.f90wrap_cputiming__get__basefnt()
    
    @basefnt.setter
    def basefnt(self, basefnt):
        _spec.f90wrap_cputiming__set__basefnt(basefnt)
    
    @property
    def tmemory(self):
        """
        Element tmemory ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 119
        
        """
        return _spec.f90wrap_cputiming__get__tmemory()
    
    @tmemory.setter
    def tmemory(self, tmemory):
        _spec.f90wrap_cputiming__set__tmemory(tmemory)
    
    @property
    def memoryt(self):
        """
        Element memoryt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 119
        
        """
        return _spec.f90wrap_cputiming__get__memoryt()
    
    @memoryt.setter
    def memoryt(self, memoryt):
        _spec.f90wrap_cputiming__set__memoryt(memoryt)
    
    @property
    def tmetrix(self):
        """
        Element tmetrix ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 120
        
        """
        return _spec.f90wrap_cputiming__get__tmetrix()
    
    @tmetrix.setter
    def tmetrix(self, tmetrix):
        _spec.f90wrap_cputiming__set__tmetrix(tmetrix)
    
    @property
    def metrixt(self):
        """
        Element metrixt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 120
        
        """
        return _spec.f90wrap_cputiming__get__metrixt()
    
    @metrixt.setter
    def metrixt(self, metrixt):
        _spec.f90wrap_cputiming__set__metrixt(metrixt)
    
    @property
    def tma00aa(self):
        """
        Element tma00aa ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 121
        
        """
        return _spec.f90wrap_cputiming__get__tma00aa()
    
    @tma00aa.setter
    def tma00aa(self, tma00aa):
        _spec.f90wrap_cputiming__set__tma00aa(tma00aa)
    
    @property
    def ma00aat(self):
        """
        Element ma00aat ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 121
        
        """
        return _spec.f90wrap_cputiming__get__ma00aat()
    
    @ma00aat.setter
    def ma00aat(self, ma00aat):
        _spec.f90wrap_cputiming__set__ma00aat(ma00aat)
    
    @property
    def tmatrix(self):
        """
        Element tmatrix ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 122
        
        """
        return _spec.f90wrap_cputiming__get__tmatrix()
    
    @tmatrix.setter
    def tmatrix(self, tmatrix):
        _spec.f90wrap_cputiming__set__tmatrix(tmatrix)
    
    @property
    def matrixt(self):
        """
        Element matrixt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 122
        
        """
        return _spec.f90wrap_cputiming__get__matrixt()
    
    @matrixt.setter
    def matrixt(self, matrixt):
        _spec.f90wrap_cputiming__set__matrixt(matrixt)
    
    @property
    def tspsmat(self):
        """
        Element tspsmat ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 123
        
        """
        return _spec.f90wrap_cputiming__get__tspsmat()
    
    @tspsmat.setter
    def tspsmat(self, tspsmat):
        _spec.f90wrap_cputiming__set__tspsmat(tspsmat)
    
    @property
    def spsmatt(self):
        """
        Element spsmatt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 123
        
        """
        return _spec.f90wrap_cputiming__get__spsmatt()
    
    @spsmatt.setter
    def spsmatt(self, spsmatt):
        _spec.f90wrap_cputiming__set__spsmatt(spsmatt)
    
    @property
    def tspsint(self):
        """
        Element tspsint ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 124
        
        """
        return _spec.f90wrap_cputiming__get__tspsint()
    
    @tspsint.setter
    def tspsint(self, tspsint):
        _spec.f90wrap_cputiming__set__tspsint(tspsint)
    
    @property
    def spsintt(self):
        """
        Element spsintt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 124
        
        """
        return _spec.f90wrap_cputiming__get__spsintt()
    
    @spsintt.setter
    def spsintt(self, spsintt):
        _spec.f90wrap_cputiming__set__spsintt(spsintt)
    
    @property
    def tmp00ac(self):
        """
        Element tmp00ac ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 125
        
        """
        return _spec.f90wrap_cputiming__get__tmp00ac()
    
    @tmp00ac.setter
    def tmp00ac(self, tmp00ac):
        _spec.f90wrap_cputiming__set__tmp00ac(tmp00ac)
    
    @property
    def mp00act(self):
        """
        Element mp00act ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 125
        
        """
        return _spec.f90wrap_cputiming__get__mp00act()
    
    @mp00act.setter
    def mp00act(self, mp00act):
        _spec.f90wrap_cputiming__set__mp00act(mp00act)
    
    @property
    def tma02aa(self):
        """
        Element tma02aa ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 126
        
        """
        return _spec.f90wrap_cputiming__get__tma02aa()
    
    @tma02aa.setter
    def tma02aa(self, tma02aa):
        _spec.f90wrap_cputiming__set__tma02aa(tma02aa)
    
    @property
    def ma02aat(self):
        """
        Element ma02aat ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 126
        
        """
        return _spec.f90wrap_cputiming__get__ma02aat()
    
    @ma02aat.setter
    def ma02aat(self, ma02aat):
        _spec.f90wrap_cputiming__set__ma02aat(ma02aat)
    
    @property
    def tpackab(self):
        """
        Element tpackab ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 127
        
        """
        return _spec.f90wrap_cputiming__get__tpackab()
    
    @tpackab.setter
    def tpackab(self, tpackab):
        _spec.f90wrap_cputiming__set__tpackab(tpackab)
    
    @property
    def packabt(self):
        """
        Element packabt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 127
        
        """
        return _spec.f90wrap_cputiming__get__packabt()
    
    @packabt.setter
    def packabt(self, packabt):
        _spec.f90wrap_cputiming__set__packabt(packabt)
    
    @property
    def ttr00ab(self):
        """
        Element ttr00ab ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 128
        
        """
        return _spec.f90wrap_cputiming__get__ttr00ab()
    
    @ttr00ab.setter
    def ttr00ab(self, ttr00ab):
        _spec.f90wrap_cputiming__set__ttr00ab(ttr00ab)
    
    @property
    def tr00abt(self):
        """
        Element tr00abt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 128
        
        """
        return _spec.f90wrap_cputiming__get__tr00abt()
    
    @tr00abt.setter
    def tr00abt(self, tr00abt):
        _spec.f90wrap_cputiming__set__tr00abt(tr00abt)
    
    @property
    def tcurent(self):
        """
        Element tcurent ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 129
        
        """
        return _spec.f90wrap_cputiming__get__tcurent()
    
    @tcurent.setter
    def tcurent(self, tcurent):
        _spec.f90wrap_cputiming__set__tcurent(tcurent)
    
    @property
    def curentt(self):
        """
        Element curentt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 129
        
        """
        return _spec.f90wrap_cputiming__get__curentt()
    
    @curentt.setter
    def curentt(self, curentt):
        _spec.f90wrap_cputiming__set__curentt(curentt)
    
    @property
    def tdf00ab(self):
        """
        Element tdf00ab ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 130
        
        """
        return _spec.f90wrap_cputiming__get__tdf00ab()
    
    @tdf00ab.setter
    def tdf00ab(self, tdf00ab):
        _spec.f90wrap_cputiming__set__tdf00ab(tdf00ab)
    
    @property
    def df00abt(self):
        """
        Element df00abt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 130
        
        """
        return _spec.f90wrap_cputiming__get__df00abt()
    
    @df00abt.setter
    def df00abt(self, df00abt):
        _spec.f90wrap_cputiming__set__df00abt(df00abt)
    
    @property
    def tlforce(self):
        """
        Element tlforce ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 131
        
        """
        return _spec.f90wrap_cputiming__get__tlforce()
    
    @tlforce.setter
    def tlforce(self, tlforce):
        _spec.f90wrap_cputiming__set__tlforce(tlforce)
    
    @property
    def lforcet(self):
        """
        Element lforcet ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 131
        
        """
        return _spec.f90wrap_cputiming__get__lforcet()
    
    @lforcet.setter
    def lforcet(self, lforcet):
        _spec.f90wrap_cputiming__set__lforcet(lforcet)
    
    @property
    def tintghs(self):
        """
        Element tintghs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 132
        
        """
        return _spec.f90wrap_cputiming__get__tintghs()
    
    @tintghs.setter
    def tintghs(self, tintghs):
        _spec.f90wrap_cputiming__set__tintghs(tintghs)
    
    @property
    def intghst(self):
        """
        Element intghst ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 132
        
        """
        return _spec.f90wrap_cputiming__get__intghst()
    
    @intghst.setter
    def intghst(self, intghst):
        _spec.f90wrap_cputiming__set__intghst(intghst)
    
    @property
    def tmtrxhs(self):
        """
        Element tmtrxhs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 133
        
        """
        return _spec.f90wrap_cputiming__get__tmtrxhs()
    
    @tmtrxhs.setter
    def tmtrxhs(self, tmtrxhs):
        _spec.f90wrap_cputiming__set__tmtrxhs(tmtrxhs)
    
    @property
    def mtrxhst(self):
        """
        Element mtrxhst ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 133
        
        """
        return _spec.f90wrap_cputiming__get__mtrxhst()
    
    @mtrxhst.setter
    def mtrxhst(self, mtrxhst):
        _spec.f90wrap_cputiming__set__mtrxhst(mtrxhst)
    
    @property
    def tlbpol(self):
        """
        Element tlbpol ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 134
        
        """
        return _spec.f90wrap_cputiming__get__tlbpol()
    
    @tlbpol.setter
    def tlbpol(self, tlbpol):
        _spec.f90wrap_cputiming__set__tlbpol(tlbpol)
    
    @property
    def lbpolt(self):
        """
        Element lbpolt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 134
        
        """
        return _spec.f90wrap_cputiming__get__lbpolt()
    
    @lbpolt.setter
    def lbpolt(self, lbpolt):
        _spec.f90wrap_cputiming__set__lbpolt(lbpolt)
    
    @property
    def tbrcast(self):
        """
        Element tbrcast ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 135
        
        """
        return _spec.f90wrap_cputiming__get__tbrcast()
    
    @tbrcast.setter
    def tbrcast(self, tbrcast):
        _spec.f90wrap_cputiming__set__tbrcast(tbrcast)
    
    @property
    def brcastt(self):
        """
        Element brcastt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 135
        
        """
        return _spec.f90wrap_cputiming__get__brcastt()
    
    @brcastt.setter
    def brcastt(self, brcastt):
        _spec.f90wrap_cputiming__set__brcastt(brcastt)
    
    @property
    def tdfp100(self):
        """
        Element tdfp100 ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 136
        
        """
        return _spec.f90wrap_cputiming__get__tdfp100()
    
    @tdfp100.setter
    def tdfp100(self, tdfp100):
        _spec.f90wrap_cputiming__set__tdfp100(tdfp100)
    
    @property
    def dfp100t(self):
        """
        Element dfp100t ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 136
        
        """
        return _spec.f90wrap_cputiming__get__dfp100t()
    
    @dfp100t.setter
    def dfp100t(self, dfp100t):
        _spec.f90wrap_cputiming__set__dfp100t(dfp100t)
    
    @property
    def tdfp200(self):
        """
        Element tdfp200 ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 137
        
        """
        return _spec.f90wrap_cputiming__get__tdfp200()
    
    @tdfp200.setter
    def tdfp200(self, tdfp200):
        _spec.f90wrap_cputiming__set__tdfp200(tdfp200)
    
    @property
    def dfp200t(self):
        """
        Element dfp200t ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 137
        
        """
        return _spec.f90wrap_cputiming__get__dfp200t()
    
    @dfp200t.setter
    def dfp200t(self, dfp200t):
        _spec.f90wrap_cputiming__set__dfp200t(dfp200t)
    
    @property
    def tdforce(self):
        """
        Element tdforce ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 138
        
        """
        return _spec.f90wrap_cputiming__get__tdforce()
    
    @tdforce.setter
    def tdforce(self, tdforce):
        _spec.f90wrap_cputiming__set__tdforce(tdforce)
    
    @property
    def dforcet(self):
        """
        Element dforcet ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 138
        
        """
        return _spec.f90wrap_cputiming__get__dforcet()
    
    @dforcet.setter
    def dforcet(self, dforcet):
        _spec.f90wrap_cputiming__set__dforcet(dforcet)
    
    @property
    def tnewton(self):
        """
        Element tnewton ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 139
        
        """
        return _spec.f90wrap_cputiming__get__tnewton()
    
    @tnewton.setter
    def tnewton(self, tnewton):
        _spec.f90wrap_cputiming__set__tnewton(tnewton)
    
    @property
    def newtont(self):
        """
        Element newtont ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 139
        
        """
        return _spec.f90wrap_cputiming__get__newtont()
    
    @newtont.setter
    def newtont(self, newtont):
        _spec.f90wrap_cputiming__set__newtont(newtont)
    
    @property
    def tcasing(self):
        """
        Element tcasing ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 140
        
        """
        return _spec.f90wrap_cputiming__get__tcasing()
    
    @tcasing.setter
    def tcasing(self, tcasing):
        _spec.f90wrap_cputiming__set__tcasing(tcasing)
    
    @property
    def casingt(self):
        """
        Element casingt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 140
        
        """
        return _spec.f90wrap_cputiming__get__casingt()
    
    @casingt.setter
    def casingt(self, casingt):
        _spec.f90wrap_cputiming__set__casingt(casingt)
    
    @property
    def tbnorml(self):
        """
        Element tbnorml ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 141
        
        """
        return _spec.f90wrap_cputiming__get__tbnorml()
    
    @tbnorml.setter
    def tbnorml(self, tbnorml):
        _spec.f90wrap_cputiming__set__tbnorml(tbnorml)
    
    @property
    def bnormlt(self):
        """
        Element bnormlt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 141
        
        """
        return _spec.f90wrap_cputiming__get__bnormlt()
    
    @bnormlt.setter
    def bnormlt(self, bnormlt):
        _spec.f90wrap_cputiming__set__bnormlt(bnormlt)
    
    @property
    def tjo00aa(self):
        """
        Element tjo00aa ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 142
        
        """
        return _spec.f90wrap_cputiming__get__tjo00aa()
    
    @tjo00aa.setter
    def tjo00aa(self, tjo00aa):
        _spec.f90wrap_cputiming__set__tjo00aa(tjo00aa)
    
    @property
    def jo00aat(self):
        """
        Element jo00aat ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 142
        
        """
        return _spec.f90wrap_cputiming__get__jo00aat()
    
    @jo00aat.setter
    def jo00aat(self, jo00aat):
        _spec.f90wrap_cputiming__set__jo00aat(jo00aat)
    
    @property
    def tpp00aa(self):
        """
        Element tpp00aa ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 143
        
        """
        return _spec.f90wrap_cputiming__get__tpp00aa()
    
    @tpp00aa.setter
    def tpp00aa(self, tpp00aa):
        _spec.f90wrap_cputiming__set__tpp00aa(tpp00aa)
    
    @property
    def pp00aat(self):
        """
        Element pp00aat ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 143
        
        """
        return _spec.f90wrap_cputiming__get__pp00aat()
    
    @pp00aat.setter
    def pp00aat(self, pp00aat):
        _spec.f90wrap_cputiming__set__pp00aat(pp00aat)
    
    @property
    def tpp00ab(self):
        """
        Element tpp00ab ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 144
        
        """
        return _spec.f90wrap_cputiming__get__tpp00ab()
    
    @tpp00ab.setter
    def tpp00ab(self, tpp00ab):
        _spec.f90wrap_cputiming__set__tpp00ab(tpp00ab)
    
    @property
    def pp00abt(self):
        """
        Element pp00abt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 144
        
        """
        return _spec.f90wrap_cputiming__get__pp00abt()
    
    @pp00abt.setter
    def pp00abt(self, pp00abt):
        _spec.f90wrap_cputiming__set__pp00abt(pp00abt)
    
    @property
    def tbfield(self):
        """
        Element tbfield ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 145
        
        """
        return _spec.f90wrap_cputiming__get__tbfield()
    
    @tbfield.setter
    def tbfield(self, tbfield):
        _spec.f90wrap_cputiming__set__tbfield(tbfield)
    
    @property
    def bfieldt(self):
        """
        Element bfieldt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 145
        
        """
        return _spec.f90wrap_cputiming__get__bfieldt()
    
    @bfieldt.setter
    def bfieldt(self, bfieldt):
        _spec.f90wrap_cputiming__set__bfieldt(bfieldt)
    
    @property
    def tstzxyz(self):
        """
        Element tstzxyz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 146
        
        """
        return _spec.f90wrap_cputiming__get__tstzxyz()
    
    @tstzxyz.setter
    def tstzxyz(self, tstzxyz):
        _spec.f90wrap_cputiming__set__tstzxyz(tstzxyz)
    
    @property
    def stzxyzt(self):
        """
        Element stzxyzt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 146
        
        """
        return _spec.f90wrap_cputiming__get__stzxyzt()
    
    @stzxyzt.setter
    def stzxyzt(self, stzxyzt):
        _spec.f90wrap_cputiming__set__stzxyzt(stzxyzt)
    
    @property
    def thesian(self):
        """
        Element thesian ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 147
        
        """
        return _spec.f90wrap_cputiming__get__thesian()
    
    @thesian.setter
    def thesian(self, thesian):
        _spec.f90wrap_cputiming__set__thesian(thesian)
    
    @property
    def hesiant(self):
        """
        Element hesiant ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 147
        
        """
        return _spec.f90wrap_cputiming__get__hesiant()
    
    @hesiant.setter
    def hesiant(self, hesiant):
        _spec.f90wrap_cputiming__set__hesiant(hesiant)
    
    @property
    def tra00aa(self):
        """
        Element tra00aa ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 148
        
        """
        return _spec.f90wrap_cputiming__get__tra00aa()
    
    @tra00aa.setter
    def tra00aa(self, tra00aa):
        _spec.f90wrap_cputiming__set__tra00aa(tra00aa)
    
    @property
    def ra00aat(self):
        """
        Element ra00aat ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 148
        
        """
        return _spec.f90wrap_cputiming__get__ra00aat()
    
    @ra00aat.setter
    def ra00aat(self, ra00aat):
        _spec.f90wrap_cputiming__set__ra00aat(ra00aat)
    
    @property
    def tnumrec(self):
        """
        Element tnumrec ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 149
        
        """
        return _spec.f90wrap_cputiming__get__tnumrec()
    
    @tnumrec.setter
    def tnumrec(self, tnumrec):
        _spec.f90wrap_cputiming__set__tnumrec(tnumrec)
    
    @property
    def numrect(self):
        """
        Element numrect ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 149
        
        """
        return _spec.f90wrap_cputiming__get__numrect()
    
    @numrect.setter
    def numrect(self, numrect):
        _spec.f90wrap_cputiming__set__numrect(numrect)
    
    @property
    def txspech(self):
        """
        Element txspech ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 150
        
        """
        return _spec.f90wrap_cputiming__get__txspech()
    
    @txspech.setter
    def txspech(self, txspech):
        _spec.f90wrap_cputiming__set__txspech(txspech)
    
    @property
    def xspecht(self):
        """
        Element xspecht ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 150
        
        """
        return _spec.f90wrap_cputiming__get__xspecht()
    
    @xspecht.setter
    def xspecht(self, xspecht):
        _spec.f90wrap_cputiming__set__xspecht(xspecht)
    
    @property
    def treadin(self):
        """
        Element treadin ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 152
        
        """
        return _spec.f90wrap_cputiming__get__treadin()
    
    @treadin.setter
    def treadin(self, treadin):
        _spec.f90wrap_cputiming__set__treadin(treadin)
    
    @property
    def twritin(self):
        """
        Element twritin ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 153
        
        """
        return _spec.f90wrap_cputiming__get__twritin()
    
    @twritin.setter
    def twritin(self, twritin):
        _spec.f90wrap_cputiming__set__twritin(twritin)
    
    @property
    def twrtend(self):
        """
        Element twrtend ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 154
        
        """
        return _spec.f90wrap_cputiming__get__twrtend()
    
    @twrtend.setter
    def twrtend(self, twrtend):
        _spec.f90wrap_cputiming__set__twrtend(twrtend)
    
    def __str__(self):
        ret = ['<cputiming>{\n']
        ret.append('    tdcuhre : ')
        ret.append(repr(self.tdcuhre))
        ret.append(',\n    dcuhret : ')
        ret.append(repr(self.dcuhret))
        ret.append(',\n    tminpack : ')
        ret.append(repr(self.tminpack))
        ret.append(',\n    minpackt : ')
        ret.append(repr(self.minpackt))
        ret.append(',\n    tiqpack : ')
        ret.append(repr(self.tiqpack))
        ret.append(',\n    iqpackt : ')
        ret.append(repr(self.iqpackt))
        ret.append(',\n    trksuite : ')
        ret.append(repr(self.trksuite))
        ret.append(',\n    rksuitet : ')
        ret.append(repr(self.rksuitet))
        ret.append(',\n    ti1mach : ')
        ret.append(repr(self.ti1mach))
        ret.append(',\n    i1macht : ')
        ret.append(repr(self.i1macht))
        ret.append(',\n    td1mach : ')
        ret.append(repr(self.td1mach))
        ret.append(',\n    d1macht : ')
        ret.append(repr(self.d1macht))
        ret.append(',\n    tilut : ')
        ret.append(repr(self.tilut))
        ret.append(',\n    ilutt : ')
        ret.append(repr(self.ilutt))
        ret.append(',\n    titers : ')
        ret.append(repr(self.titers))
        ret.append(',\n    iterst : ')
        ret.append(repr(self.iterst))
        ret.append(',\n    tinputlist : ')
        ret.append(repr(self.tinputlist))
        ret.append(',\n    inputlistt : ')
        ret.append(repr(self.inputlistt))
        ret.append(',\n    tglobal : ')
        ret.append(repr(self.tglobal))
        ret.append(',\n    globalt : ')
        ret.append(repr(self.globalt))
        ret.append(',\n    tsphdf5 : ')
        ret.append(repr(self.tsphdf5))
        ret.append(',\n    sphdf5t : ')
        ret.append(repr(self.sphdf5t))
        ret.append(',\n    tpreset : ')
        ret.append(repr(self.tpreset))
        ret.append(',\n    presett : ')
        ret.append(repr(self.presett))
        ret.append(',\n    tmanual : ')
        ret.append(repr(self.tmanual))
        ret.append(',\n    manualt : ')
        ret.append(repr(self.manualt))
        ret.append(',\n    trzaxis : ')
        ret.append(repr(self.trzaxis))
        ret.append(',\n    rzaxist : ')
        ret.append(repr(self.rzaxist))
        ret.append(',\n    tpackxi : ')
        ret.append(repr(self.tpackxi))
        ret.append(',\n    packxit : ')
        ret.append(repr(self.packxit))
        ret.append(',\n    tvolume : ')
        ret.append(repr(self.tvolume))
        ret.append(',\n    volumet : ')
        ret.append(repr(self.volumet))
        ret.append(',\n    tcoords : ')
        ret.append(repr(self.tcoords))
        ret.append(',\n    coordst : ')
        ret.append(repr(self.coordst))
        ret.append(',\n    tbasefn : ')
        ret.append(repr(self.tbasefn))
        ret.append(',\n    basefnt : ')
        ret.append(repr(self.basefnt))
        ret.append(',\n    tmemory : ')
        ret.append(repr(self.tmemory))
        ret.append(',\n    memoryt : ')
        ret.append(repr(self.memoryt))
        ret.append(',\n    tmetrix : ')
        ret.append(repr(self.tmetrix))
        ret.append(',\n    metrixt : ')
        ret.append(repr(self.metrixt))
        ret.append(',\n    tma00aa : ')
        ret.append(repr(self.tma00aa))
        ret.append(',\n    ma00aat : ')
        ret.append(repr(self.ma00aat))
        ret.append(',\n    tmatrix : ')
        ret.append(repr(self.tmatrix))
        ret.append(',\n    matrixt : ')
        ret.append(repr(self.matrixt))
        ret.append(',\n    tspsmat : ')
        ret.append(repr(self.tspsmat))
        ret.append(',\n    spsmatt : ')
        ret.append(repr(self.spsmatt))
        ret.append(',\n    tspsint : ')
        ret.append(repr(self.tspsint))
        ret.append(',\n    spsintt : ')
        ret.append(repr(self.spsintt))
        ret.append(',\n    tmp00ac : ')
        ret.append(repr(self.tmp00ac))
        ret.append(',\n    mp00act : ')
        ret.append(repr(self.mp00act))
        ret.append(',\n    tma02aa : ')
        ret.append(repr(self.tma02aa))
        ret.append(',\n    ma02aat : ')
        ret.append(repr(self.ma02aat))
        ret.append(',\n    tpackab : ')
        ret.append(repr(self.tpackab))
        ret.append(',\n    packabt : ')
        ret.append(repr(self.packabt))
        ret.append(',\n    ttr00ab : ')
        ret.append(repr(self.ttr00ab))
        ret.append(',\n    tr00abt : ')
        ret.append(repr(self.tr00abt))
        ret.append(',\n    tcurent : ')
        ret.append(repr(self.tcurent))
        ret.append(',\n    curentt : ')
        ret.append(repr(self.curentt))
        ret.append(',\n    tdf00ab : ')
        ret.append(repr(self.tdf00ab))
        ret.append(',\n    df00abt : ')
        ret.append(repr(self.df00abt))
        ret.append(',\n    tlforce : ')
        ret.append(repr(self.tlforce))
        ret.append(',\n    lforcet : ')
        ret.append(repr(self.lforcet))
        ret.append(',\n    tintghs : ')
        ret.append(repr(self.tintghs))
        ret.append(',\n    intghst : ')
        ret.append(repr(self.intghst))
        ret.append(',\n    tmtrxhs : ')
        ret.append(repr(self.tmtrxhs))
        ret.append(',\n    mtrxhst : ')
        ret.append(repr(self.mtrxhst))
        ret.append(',\n    tlbpol : ')
        ret.append(repr(self.tlbpol))
        ret.append(',\n    lbpolt : ')
        ret.append(repr(self.lbpolt))
        ret.append(',\n    tbrcast : ')
        ret.append(repr(self.tbrcast))
        ret.append(',\n    brcastt : ')
        ret.append(repr(self.brcastt))
        ret.append(',\n    tdfp100 : ')
        ret.append(repr(self.tdfp100))
        ret.append(',\n    dfp100t : ')
        ret.append(repr(self.dfp100t))
        ret.append(',\n    tdfp200 : ')
        ret.append(repr(self.tdfp200))
        ret.append(',\n    dfp200t : ')
        ret.append(repr(self.dfp200t))
        ret.append(',\n    tdforce : ')
        ret.append(repr(self.tdforce))
        ret.append(',\n    dforcet : ')
        ret.append(repr(self.dforcet))
        ret.append(',\n    tnewton : ')
        ret.append(repr(self.tnewton))
        ret.append(',\n    newtont : ')
        ret.append(repr(self.newtont))
        ret.append(',\n    tcasing : ')
        ret.append(repr(self.tcasing))
        ret.append(',\n    casingt : ')
        ret.append(repr(self.casingt))
        ret.append(',\n    tbnorml : ')
        ret.append(repr(self.tbnorml))
        ret.append(',\n    bnormlt : ')
        ret.append(repr(self.bnormlt))
        ret.append(',\n    tjo00aa : ')
        ret.append(repr(self.tjo00aa))
        ret.append(',\n    jo00aat : ')
        ret.append(repr(self.jo00aat))
        ret.append(',\n    tpp00aa : ')
        ret.append(repr(self.tpp00aa))
        ret.append(',\n    pp00aat : ')
        ret.append(repr(self.pp00aat))
        ret.append(',\n    tpp00ab : ')
        ret.append(repr(self.tpp00ab))
        ret.append(',\n    pp00abt : ')
        ret.append(repr(self.pp00abt))
        ret.append(',\n    tbfield : ')
        ret.append(repr(self.tbfield))
        ret.append(',\n    bfieldt : ')
        ret.append(repr(self.bfieldt))
        ret.append(',\n    tstzxyz : ')
        ret.append(repr(self.tstzxyz))
        ret.append(',\n    stzxyzt : ')
        ret.append(repr(self.stzxyzt))
        ret.append(',\n    thesian : ')
        ret.append(repr(self.thesian))
        ret.append(',\n    hesiant : ')
        ret.append(repr(self.hesiant))
        ret.append(',\n    tra00aa : ')
        ret.append(repr(self.tra00aa))
        ret.append(',\n    ra00aat : ')
        ret.append(repr(self.ra00aat))
        ret.append(',\n    tnumrec : ')
        ret.append(repr(self.tnumrec))
        ret.append(',\n    numrect : ')
        ret.append(repr(self.numrect))
        ret.append(',\n    txspech : ')
        ret.append(repr(self.txspech))
        ret.append(',\n    xspecht : ')
        ret.append(repr(self.xspecht))
        ret.append(',\n    treadin : ')
        ret.append(repr(self.treadin))
        ret.append(',\n    twritin : ')
        ret.append(repr(self.twritin))
        ret.append(',\n    twrtend : ')
        ret.append(repr(self.twrtend))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

cputiming = Cputiming()

class Typedefns(f90wrap.runtime.FortranModule):
    """
    Module typedefns
    
    
    Defined at global.fpp lines 157-173
    
    """
    @f90wrap.runtime.register_class("spec.subgrid")
    class subgrid(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=subgrid)
        
        
        Defined at global.fpp lines 158-160
        
        """
        def __init__(self, handle=None):
            """
            self = Subgrid()
            
            
            Defined at global.fpp lines 158-160
            
            
            Returns
            -------
            this : Subgrid
            	Object to be constructed
            
            
            Automatically generated constructor for subgrid
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _spec.f90wrap_subgrid_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Subgrid
            
            
            Defined at global.fpp lines 158-160
            
            Parameters
            ----------
            this : Subgrid
            	Object to be destructed
            
            
            Automatically generated destructor for subgrid
            """
            if self._alloc:
                _spec.f90wrap_subgrid_finalise(this=self._handle)
        
        @property
        def s(self):
            """
            Element s ftype=real(8) pytype=float
            
            
            Defined at global.fpp line 159
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_subgrid__array__s(self._handle)
            if array_handle in self._arrays:
                s = self._arrays[array_handle]
            else:
                s = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_subgrid__array__s)
                self._arrays[array_handle] = s
            return s
        
        @s.setter
        def s(self, s):
            self.s[...] = s
        
        @property
        def i(self):
            """
            Element i ftype=integer pytype=int
            
            
            Defined at global.fpp line 160
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_subgrid__array__i(self._handle)
            if array_handle in self._arrays:
                i = self._arrays[array_handle]
            else:
                i = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_subgrid__array__i)
                self._arrays[array_handle] = i
            return i
        
        @i.setter
        def i(self, i):
            self.i[...] = i
        
        def __str__(self):
            ret = ['<subgrid>{\n']
            ret.append('    s : ')
            ret.append(repr(self.s))
            ret.append(',\n    i : ')
            ret.append(repr(self.i))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("spec.MatrixLU")
    class MatrixLU(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=matrixlu)
        
        
        Defined at global.fpp lines 162-164
        
        """
        def __init__(self, handle=None):
            """
            self = Matrixlu()
            
            
            Defined at global.fpp lines 162-164
            
            
            Returns
            -------
            this : Matrixlu
            	Object to be constructed
            
            
            Automatically generated constructor for matrixlu
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _spec.f90wrap_matrixlu_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Matrixlu
            
            
            Defined at global.fpp lines 162-164
            
            Parameters
            ----------
            this : Matrixlu
            	Object to be destructed
            
            
            Automatically generated destructor for matrixlu
            """
            if self._alloc:
                _spec.f90wrap_matrixlu_finalise(this=self._handle)
        
        @property
        def mat(self):
            """
            Element mat ftype=real(8) pytype=float
            
            
            Defined at global.fpp line 163
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_matrixlu__array__mat(self._handle)
            if array_handle in self._arrays:
                mat = self._arrays[array_handle]
            else:
                mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_matrixlu__array__mat)
                self._arrays[array_handle] = mat
            return mat
        
        @mat.setter
        def mat(self, mat):
            self.mat[...] = mat
        
        @property
        def ipivot(self):
            """
            Element ipivot ftype=integer pytype=int
            
            
            Defined at global.fpp line 164
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_matrixlu__array__ipivot(self._handle)
            if array_handle in self._arrays:
                ipivot = self._arrays[array_handle]
            else:
                ipivot = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_matrixlu__array__ipivot)
                self._arrays[array_handle] = ipivot
            return ipivot
        
        @ipivot.setter
        def ipivot(self, ipivot):
            self.ipivot[...] = ipivot
        
        def __str__(self):
            ret = ['<matrixlu>{\n']
            ret.append('    mat : ')
            ret.append(repr(self.mat))
            ret.append(',\n    ipivot : ')
            ret.append(repr(self.ipivot))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("spec.derivative")
    class derivative(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=derivative)
        
        
        Defined at global.fpp lines 166-172
        
        """
        def __init__(self, handle=None):
            """
            self = Derivative()
            
            
            Defined at global.fpp lines 166-172
            
            
            Returns
            -------
            this : Derivative
            	Object to be constructed
            
            
            Automatically generated constructor for derivative
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _spec.f90wrap_derivative_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Derivative
            
            
            Defined at global.fpp lines 166-172
            
            Parameters
            ----------
            this : Derivative
            	Object to be destructed
            
            
            Automatically generated destructor for derivative
            """
            if self._alloc:
                _spec.f90wrap_derivative_finalise(this=self._handle)
        
        @property
        def l(self):
            """
            Element l ftype=logical pytype=bool
            
            
            Defined at global.fpp line 167
            
            """
            return _spec.f90wrap_derivative__get__l(self._handle)
        
        @l.setter
        def l(self, l):
            _spec.f90wrap_derivative__set__l(self._handle, l)
        
        @property
        def vol(self):
            """
            Element vol ftype=integer  pytype=int
            
            
            Defined at global.fpp line 168
            
            """
            return _spec.f90wrap_derivative__get__vol(self._handle)
        
        @vol.setter
        def vol(self, vol):
            _spec.f90wrap_derivative__set__vol(self._handle, vol)
        
        @property
        def innout(self):
            """
            Element innout ftype=integer  pytype=int
            
            
            Defined at global.fpp line 169
            
            """
            return _spec.f90wrap_derivative__get__innout(self._handle)
        
        @innout.setter
        def innout(self, innout):
            _spec.f90wrap_derivative__set__innout(self._handle, innout)
        
        @property
        def ii(self):
            """
            Element ii ftype=integer  pytype=int
            
            
            Defined at global.fpp line 170
            
            """
            return _spec.f90wrap_derivative__get__ii(self._handle)
        
        @ii.setter
        def ii(self, ii):
            _spec.f90wrap_derivative__set__ii(self._handle, ii)
        
        @property
        def irz(self):
            """
            Element irz ftype=integer  pytype=int
            
            
            Defined at global.fpp line 171
            
            """
            return _spec.f90wrap_derivative__get__irz(self._handle)
        
        @irz.setter
        def irz(self, irz):
            _spec.f90wrap_derivative__set__irz(self._handle, irz)
        
        @property
        def issym(self):
            """
            Element issym ftype=integer  pytype=int
            
            
            Defined at global.fpp line 172
            
            """
            return _spec.f90wrap_derivative__get__issym(self._handle)
        
        @issym.setter
        def issym(self, issym):
            _spec.f90wrap_derivative__set__issym(self._handle, issym)
        
        def __str__(self):
            ret = ['<derivative>{\n']
            ret.append('    l : ')
            ret.append(repr(self.l))
            ret.append(',\n    vol : ')
            ret.append(repr(self.vol))
            ret.append(',\n    innout : ')
            ret.append(repr(self.innout))
            ret.append(',\n    ii : ')
            ret.append(repr(self.ii))
            ret.append(',\n    irz : ')
            ret.append(repr(self.irz))
            ret.append(',\n    issym : ')
            ret.append(repr(self.issym))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

typedefns = Typedefns()

class Allglobal(f90wrap.runtime.FortranModule):
    """
    Module allglobal
    
    
    Defined at global.fpp lines 176-2271
    
    """
    @staticmethod
    def build_vector_potential(lvol, iocons, aderiv, tderiv):
        """
        build_vector_potential(lvol, iocons, aderiv, tderiv)
        
        
        Defined at global.fpp lines 536-589
        
        Parameters
        ----------
        lvol : int
        iocons : int
        aderiv : int
        tderiv : int
        
        """
        _spec.f90wrap_build_vector_potential(lvol=lvol, iocons=iocons, aderiv=aderiv, \
            tderiv=tderiv)
    
    @staticmethod
    def set_mpi_comm(comm):
        """
        set_mpi_comm(comm)
        
        
        Defined at global.fpp lines 595-607
        
        Parameters
        ----------
        comm : int
        
        """
        _spec.f90wrap_set_mpi_comm(comm=comm)
    
    @staticmethod
    def read_inputlists_from_file():
        """
        read_inputlists_from_file()
        
        
        Defined at global.fpp lines 610-733
        
        
        """
        _spec.f90wrap_read_inputlists_from_file()
    
    @staticmethod
    def check_inputs():
        """
        check_inputs()
        
        
        Defined at global.fpp lines 737-1126
        
        
        """
        _spec.f90wrap_check_inputs()
    
    @staticmethod
    def broadcast_inputs():
        """
        broadcast_inputs()
        
        
        Defined at global.fpp lines 1130-1927
        
        
        """
        _spec.f90wrap_broadcast_inputs()
    
    @staticmethod
    def wrtend():
        """
        wrtend()
        
        
        Defined at global.fpp lines 1931-2235
        
        
        """
        _spec.f90wrap_wrtend()
    
    @staticmethod
    def ismyvolume(vvol):
        """
        ismyvolume(vvol)
        
        
        Defined at global.fpp lines 2238-2255
        
        Parameters
        ----------
        vvol : int
        
        """
        _spec.f90wrap_ismyvolume(vvol=vvol)
    
    @staticmethod
    def whichcpuid(vvol, cpu_id):
        """
        whichcpuid(vvol, cpu_id)
        
        
        Defined at global.fpp lines 2258-2269
        
        Parameters
        ----------
        vvol : int
        cpu_id : int
        
        """
        _spec.f90wrap_whichcpuid(vvol=vvol, cpu_id=cpu_id)
    
    @property
    def myid(self):
        """
        Element myid ftype=integer               pytype=int
        
        
        Defined at global.fpp line 181
        
        """
        return _spec.f90wrap_allglobal__get__myid()
    
    @myid.setter
    def myid(self, myid):
        _spec.f90wrap_allglobal__set__myid(myid)
    
    @property
    def ncpu(self):
        """
        Element ncpu ftype=integer               pytype=int
        
        
        Defined at global.fpp line 181
        
        """
        return _spec.f90wrap_allglobal__get__ncpu()
    
    @ncpu.setter
    def ncpu(self, ncpu):
        _spec.f90wrap_allglobal__set__ncpu(ncpu)
    
    @property
    def mpi_comm_spec(self):
        """
        Element mpi_comm_spec ftype=integer               pytype=int
        
        
        Defined at global.fpp line 181
        
        """
        return _spec.f90wrap_allglobal__get__mpi_comm_spec()
    
    @mpi_comm_spec.setter
    def mpi_comm_spec(self, mpi_comm_spec):
        _spec.f90wrap_allglobal__set__mpi_comm_spec(mpi_comm_spec)
    
    @property
    def ismyvolumevalue(self):
        """
        Element ismyvolumevalue ftype=integer               pytype=int
        
        
        Defined at global.fpp line 182
        
        """
        return _spec.f90wrap_allglobal__get__ismyvolumevalue()
    
    @ismyvolumevalue.setter
    def ismyvolumevalue(self, ismyvolumevalue):
        _spec.f90wrap_allglobal__set__ismyvolumevalue(ismyvolumevalue)
    
    @property
    def cpus(self):
        """
        Element cpus ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 183
        
        """
        return _spec.f90wrap_allglobal__get__cpus()
    
    @cpus.setter
    def cpus(self, cpus):
        _spec.f90wrap_allglobal__set__cpus(cpus)
    
    @property
    def skip_write(self):
        """
        Element skip_write ftype=logical pytype=bool
        
        
        Defined at global.fpp line 184
        
        """
        return _spec.f90wrap_allglobal__get__skip_write()
    
    @skip_write.setter
    def skip_write(self, skip_write):
        _spec.f90wrap_allglobal__set__skip_write(skip_write)
    
    @property
    def pi2nfp(self):
        """
        Element pi2nfp ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 185
        
        """
        return _spec.f90wrap_allglobal__get__pi2nfp()
    
    @pi2nfp.setter
    def pi2nfp(self, pi2nfp):
        _spec.f90wrap_allglobal__set__pi2nfp(pi2nfp)
    
    @property
    def pi2pi2nfp(self):
        """
        Element pi2pi2nfp ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 186
        
        """
        return _spec.f90wrap_allglobal__get__pi2pi2nfp()
    
    @pi2pi2nfp.setter
    def pi2pi2nfp(self, pi2pi2nfp):
        _spec.f90wrap_allglobal__set__pi2pi2nfp(pi2pi2nfp)
    
    @property
    def pi2pi2nfphalf(self):
        """
        Element pi2pi2nfphalf ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 187
        
        """
        return _spec.f90wrap_allglobal__get__pi2pi2nfphalf()
    
    @pi2pi2nfphalf.setter
    def pi2pi2nfphalf(self, pi2pi2nfphalf):
        _spec.f90wrap_allglobal__set__pi2pi2nfphalf(pi2pi2nfphalf)
    
    @property
    def pi2pi2nfpquart(self):
        """
        Element pi2pi2nfpquart ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 188
        
        """
        return _spec.f90wrap_allglobal__get__pi2pi2nfpquart()
    
    @pi2pi2nfpquart.setter
    def pi2pi2nfpquart(self, pi2pi2nfpquart):
        _spec.f90wrap_allglobal__set__pi2pi2nfpquart(pi2pi2nfpquart)
    
    @property
    def ext(self):
        """
        Element ext ftype=character(len=100) pytype=str
        
        
        Defined at global.fpp line 190
        
        """
        return _spec.f90wrap_allglobal__get__ext()
    
    @ext.setter
    def ext(self, ext):
        _spec.f90wrap_allglobal__set__ext(ext)
    
    @property
    def forceerr(self):
        """
        Element forceerr ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 191
        
        """
        return _spec.f90wrap_allglobal__get__forceerr()
    
    @forceerr.setter
    def forceerr(self, forceerr):
        _spec.f90wrap_allglobal__set__forceerr(forceerr)
    
    @property
    def energy(self):
        """
        Element energy ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 191
        
        """
        return _spec.f90wrap_allglobal__get__energy()
    
    @energy.setter
    def energy(self, energy):
        _spec.f90wrap_allglobal__set__energy(energy)
    
    @property
    def ipdt(self):
        """
        Element ipdt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 192
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ipdt(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ipdt = self._arrays[array_handle]
        else:
            ipdt = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ipdt)
            self._arrays[array_handle] = ipdt
        return ipdt
    
    @ipdt.setter
    def ipdt(self, ipdt):
        self.ipdt[...] = ipdt
    
    @property
    def ipdtdpf(self):
        """
        Element ipdtdpf ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 192
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ipdtdpf(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ipdtdpf = self._arrays[array_handle]
        else:
            ipdtdpf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ipdtdpf)
            self._arrays[array_handle] = ipdtdpf
        return ipdtdpf
    
    @ipdtdpf.setter
    def ipdtdpf(self, ipdtdpf):
        self.ipdtdpf[...] = ipdtdpf
    
    @property
    def mvol(self):
        """
        Element mvol ftype=integer               pytype=int
        
        
        Defined at global.fpp line 193
        
        """
        return _spec.f90wrap_allglobal__get__mvol()
    
    @mvol.setter
    def mvol(self, mvol):
        _spec.f90wrap_allglobal__set__mvol(mvol)
    
    @property
    def yesstellsym(self):
        """
        Element yesstellsym ftype=logical pytype=bool
        
        
        Defined at global.fpp line 194
        
        """
        return _spec.f90wrap_allglobal__get__yesstellsym()
    
    @yesstellsym.setter
    def yesstellsym(self, yesstellsym):
        _spec.f90wrap_allglobal__set__yesstellsym(yesstellsym)
    
    @property
    def notstellsym(self):
        """
        Element notstellsym ftype=logical pytype=bool
        
        
        Defined at global.fpp line 194
        
        """
        return _spec.f90wrap_allglobal__get__notstellsym()
    
    @notstellsym.setter
    def notstellsym(self, notstellsym):
        _spec.f90wrap_allglobal__set__notstellsym(notstellsym)
    
    @property
    def yesmatrixfree(self):
        """
        Element yesmatrixfree ftype=logical pytype=bool
        
        
        Defined at global.fpp line 195
        
        """
        return _spec.f90wrap_allglobal__get__yesmatrixfree()
    
    @yesmatrixfree.setter
    def yesmatrixfree(self, yesmatrixfree):
        _spec.f90wrap_allglobal__set__yesmatrixfree(yesmatrixfree)
    
    @property
    def notmatrixfree(self):
        """
        Element notmatrixfree ftype=logical pytype=bool
        
        
        Defined at global.fpp line 195
        
        """
        return _spec.f90wrap_allglobal__get__notmatrixfree()
    
    @notmatrixfree.setter
    def notmatrixfree(self, notmatrixfree):
        _spec.f90wrap_allglobal__set__notmatrixfree(notmatrixfree)
    
    @property
    def cheby(self):
        """
        Element cheby ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__cheby(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cheby = self._arrays[array_handle]
        else:
            cheby = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__cheby)
            self._arrays[array_handle] = cheby
        return cheby
    
    @cheby.setter
    def cheby(self, cheby):
        self.cheby[...] = cheby
    
    @property
    def zernike(self):
        """
        Element zernike ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__zernike(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zernike = self._arrays[array_handle]
        else:
            zernike = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__zernike)
            self._arrays[array_handle] = zernike
        return zernike
    
    @zernike.setter
    def zernike(self, zernike):
        self.zernike[...] = zernike
    
    @property
    def tt(self):
        """
        Element tt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 197
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tt(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tt = self._arrays[array_handle]
        else:
            tt = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tt)
            self._arrays[array_handle] = tt
        return tt
    
    @tt.setter
    def tt(self, tt):
        self.tt[...] = tt
    
    @property
    def rtt(self):
        """
        Element rtt ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 197
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__rtt(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rtt = self._arrays[array_handle]
        else:
            rtt = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__rtt)
            self._arrays[array_handle] = rtt
        return rtt
    
    @rtt.setter
    def rtt(self, rtt):
        self.rtt[...] = rtt
    
    @property
    def rtm(self):
        """
        Element rtm ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 198
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__rtm(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rtm = self._arrays[array_handle]
        else:
            rtm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__rtm)
            self._arrays[array_handle] = rtm
        return rtm
    
    @rtm.setter
    def rtm(self, rtm):
        self.rtm[...] = rtm
    
    @property
    def zernikedof(self):
        """
        Element zernikedof ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 199
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__zernikedof(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zernikedof = self._arrays[array_handle]
        else:
            zernikedof = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__zernikedof)
            self._arrays[array_handle] = zernikedof
        return zernikedof
    
    @zernikedof.setter
    def zernikedof(self, zernikedof):
        self.zernikedof[...] = zernikedof
    
    @property
    def imagneticok(self):
        """
        Element imagneticok ftype=logical pytype=bool
        
        
        Defined at global.fpp line 200
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__imagneticok(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            imagneticok = self._arrays[array_handle]
        else:
            imagneticok = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__imagneticok)
            self._arrays[array_handle] = imagneticok
        return imagneticok
    
    @imagneticok.setter
    def imagneticok(self, imagneticok):
        self.imagneticok[...] = imagneticok
    
    @property
    def iconstraintok(self):
        """
        Element iconstraintok ftype=logical pytype=bool
        
        
        Defined at global.fpp line 201
        
        """
        return _spec.f90wrap_allglobal__get__iconstraintok()
    
    @iconstraintok.setter
    def iconstraintok(self, iconstraintok):
        _spec.f90wrap_allglobal__set__iconstraintok(iconstraintok)
    
    @property
    def beltramierror(self):
        """
        Element beltramierror ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 202
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__beltramierror(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            beltramierror = self._arrays[array_handle]
        else:
            beltramierror = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__beltramierror)
            self._arrays[array_handle] = beltramierror
        return beltramierror
    
    @beltramierror.setter
    def beltramierror(self, beltramierror):
        self.beltramierror[...] = beltramierror
    
    @property
    def mn(self):
        """
        Element mn ftype=integer               pytype=int
        
        
        Defined at global.fpp line 207
        
        """
        return _spec.f90wrap_allglobal__get__mn()
    
    @mn.setter
    def mn(self, mn):
        _spec.f90wrap_allglobal__set__mn(mn)
    
    @property
    def im(self):
        """
        Element im ftype=integer pytype=int
        
        
        Defined at global.fpp line 208
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__im(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            im = self._arrays[array_handle]
        else:
            im = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__im)
            self._arrays[array_handle] = im
        return im
    
    @im.setter
    def im(self, im):
        self.im[...] = im
    
    @property
    def in_(self):
        """
        Element in_ ftype=integer pytype=int
        
        
        Defined at global.fpp line 208
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__in_(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            in_ = self._arrays[array_handle]
        else:
            in_ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__in_)
            self._arrays[array_handle] = in_
        return in_
    
    @in_.setter
    def in_(self, in_):
        self.in_[...] = in_
    
    @property
    def halfmm(self):
        """
        Element halfmm ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 209
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__halfmm(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            halfmm = self._arrays[array_handle]
        else:
            halfmm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__halfmm)
            self._arrays[array_handle] = halfmm
        return halfmm
    
    @halfmm.setter
    def halfmm(self, halfmm):
        self.halfmm[...] = halfmm
    
    @property
    def regumm(self):
        """
        Element regumm ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 209
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__regumm(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            regumm = self._arrays[array_handle]
        else:
            regumm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__regumm)
            self._arrays[array_handle] = regumm
        return regumm
    
    @regumm.setter
    def regumm(self, regumm):
        self.regumm[...] = regumm
    
    @property
    def rscale(self):
        """
        Element rscale ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 210
        
        """
        return _spec.f90wrap_allglobal__get__rscale()
    
    @rscale.setter
    def rscale(self, rscale):
        _spec.f90wrap_allglobal__set__rscale(rscale)
    
    @property
    def psifactor(self):
        """
        Element psifactor ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 211
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__psifactor(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psifactor = self._arrays[array_handle]
        else:
            psifactor = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__psifactor)
            self._arrays[array_handle] = psifactor
        return psifactor
    
    @psifactor.setter
    def psifactor(self, psifactor):
        self.psifactor[...] = psifactor
    
    @property
    def inifactor(self):
        """
        Element inifactor ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 211
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__inifactor(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            inifactor = self._arrays[array_handle]
        else:
            inifactor = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__inifactor)
            self._arrays[array_handle] = inifactor
        return inifactor
    
    @inifactor.setter
    def inifactor(self, inifactor):
        self.inifactor[...] = inifactor
    
    @property
    def bbweight(self):
        """
        Element bbweight ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 212
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bbweight(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bbweight = self._arrays[array_handle]
        else:
            bbweight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bbweight)
            self._arrays[array_handle] = bbweight
        return bbweight
    
    @bbweight.setter
    def bbweight(self, bbweight):
        self.bbweight[...] = bbweight
    
    @property
    def mmpp(self):
        """
        Element mmpp ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 213
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__mmpp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            mmpp = self._arrays[array_handle]
        else:
            mmpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__mmpp)
            self._arrays[array_handle] = mmpp
        return mmpp
    
    @mmpp.setter
    def mmpp(self, mmpp):
        self.mmpp[...] = mmpp
    
    @property
    def mne(self):
        """
        Element mne ftype=integer               pytype=int
        
        
        Defined at global.fpp line 217
        
        """
        return _spec.f90wrap_allglobal__get__mne()
    
    @mne.setter
    def mne(self, mne):
        _spec.f90wrap_allglobal__set__mne(mne)
    
    @property
    def ime(self):
        """
        Element ime ftype=integer pytype=int
        
        
        Defined at global.fpp line 218
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ime(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ime = self._arrays[array_handle]
        else:
            ime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ime)
            self._arrays[array_handle] = ime
        return ime
    
    @ime.setter
    def ime(self, ime):
        self.ime[...] = ime
    
    @property
    def ine(self):
        """
        Element ine ftype=integer pytype=int
        
        
        Defined at global.fpp line 218
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ine(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ine = self._arrays[array_handle]
        else:
            ine = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ine)
            self._arrays[array_handle] = ine
        return ine
    
    @ine.setter
    def ine(self, ine):
        self.ine[...] = ine
    
    @property
    def mns(self):
        """
        Element mns ftype=integer               pytype=int
        
        
        Defined at global.fpp line 223
        
        """
        return _spec.f90wrap_allglobal__get__mns()
    
    @mns.setter
    def mns(self, mns):
        _spec.f90wrap_allglobal__set__mns(mns)
    
    @property
    def ims(self):
        """
        Element ims ftype=integer pytype=int
        
        
        Defined at global.fpp line 224
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ims(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ims = self._arrays[array_handle]
        else:
            ims = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ims)
            self._arrays[array_handle] = ims
        return ims
    
    @ims.setter
    def ims(self, ims):
        self.ims[...] = ims
    
    @property
    def ins(self):
        """
        Element ins ftype=integer pytype=int
        
        
        Defined at global.fpp line 224
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ins(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ins = self._arrays[array_handle]
        else:
            ins = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ins)
            self._arrays[array_handle] = ins
        return ins
    
    @ins.setter
    def ins(self, ins):
        self.ins[...] = ins
    
    @property
    def lmpol(self):
        """
        Element lmpol ftype=integer               pytype=int
        
        
        Defined at global.fpp line 225
        
        """
        return _spec.f90wrap_allglobal__get__lmpol()
    
    @lmpol.setter
    def lmpol(self, lmpol):
        _spec.f90wrap_allglobal__set__lmpol(lmpol)
    
    @property
    def lntor(self):
        """
        Element lntor ftype=integer               pytype=int
        
        
        Defined at global.fpp line 225
        
        """
        return _spec.f90wrap_allglobal__get__lntor()
    
    @lntor.setter
    def lntor(self, lntor):
        _spec.f90wrap_allglobal__set__lntor(lntor)
    
    @property
    def smpol(self):
        """
        Element smpol ftype=integer               pytype=int
        
        
        Defined at global.fpp line 225
        
        """
        return _spec.f90wrap_allglobal__get__smpol()
    
    @smpol.setter
    def smpol(self, smpol):
        _spec.f90wrap_allglobal__set__smpol(smpol)
    
    @property
    def sntor(self):
        """
        Element sntor ftype=integer               pytype=int
        
        
        Defined at global.fpp line 225
        
        """
        return _spec.f90wrap_allglobal__get__sntor()
    
    @sntor.setter
    def sntor(self, sntor):
        _spec.f90wrap_allglobal__set__sntor(sntor)
    
    @property
    def xoffset(self):
        """
        Element xoffset ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 227
        
        """
        return _spec.f90wrap_allglobal__get__xoffset()
    
    @xoffset.setter
    def xoffset(self, xoffset):
        _spec.f90wrap_allglobal__set__xoffset(xoffset)
    
    @property
    def irbc(self):
        """
        Element irbc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 234
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__irbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            irbc = self._arrays[array_handle]
        else:
            irbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__irbc)
            self._arrays[array_handle] = irbc
        return irbc
    
    @irbc.setter
    def irbc(self, irbc):
        self.irbc[...] = irbc
    
    @property
    def izbs(self):
        """
        Element izbs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 234
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__izbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            izbs = self._arrays[array_handle]
        else:
            izbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__izbs)
            self._arrays[array_handle] = izbs
        return izbs
    
    @izbs.setter
    def izbs(self, izbs):
        self.izbs[...] = izbs
    
    @property
    def irbs(self):
        """
        Element irbs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 235
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__irbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            irbs = self._arrays[array_handle]
        else:
            irbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__irbs)
            self._arrays[array_handle] = irbs
        return irbs
    
    @irbs.setter
    def irbs(self, irbs):
        self.irbs[...] = irbs
    
    @property
    def izbc(self):
        """
        Element izbc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 235
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__izbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            izbc = self._arrays[array_handle]
        else:
            izbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__izbc)
            self._arrays[array_handle] = izbc
        return izbc
    
    @izbc.setter
    def izbc(self, izbc):
        self.izbc[...] = izbc
    
    @property
    def drbc(self):
        """
        Element drbc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 236
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__drbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            drbc = self._arrays[array_handle]
        else:
            drbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__drbc)
            self._arrays[array_handle] = drbc
        return drbc
    
    @drbc.setter
    def drbc(self, drbc):
        self.drbc[...] = drbc
    
    @property
    def dzbs(self):
        """
        Element dzbs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 236
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzbs = self._arrays[array_handle]
        else:
            dzbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzbs)
            self._arrays[array_handle] = dzbs
        return dzbs
    
    @dzbs.setter
    def dzbs(self, dzbs):
        self.dzbs[...] = dzbs
    
    @property
    def drbs(self):
        """
        Element drbs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 237
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__drbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            drbs = self._arrays[array_handle]
        else:
            drbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__drbs)
            self._arrays[array_handle] = drbs
        return drbs
    
    @drbs.setter
    def drbs(self, drbs):
        self.drbs[...] = drbs
    
    @property
    def dzbc(self):
        """
        Element dzbc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 237
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzbc = self._arrays[array_handle]
        else:
            dzbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzbc)
            self._arrays[array_handle] = dzbc
        return dzbc
    
    @dzbc.setter
    def dzbc(self, dzbc):
        self.dzbc[...] = dzbc
    
    @property
    def irij(self):
        """
        Element irij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 238
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__irij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            irij = self._arrays[array_handle]
        else:
            irij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__irij)
            self._arrays[array_handle] = irij
        return irij
    
    @irij.setter
    def irij(self, irij):
        self.irij[...] = irij
    
    @property
    def izij(self):
        """
        Element izij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 238
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__izij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            izij = self._arrays[array_handle]
        else:
            izij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__izij)
            self._arrays[array_handle] = izij
        return izij
    
    @izij.setter
    def izij(self, izij):
        self.izij[...] = izij
    
    @property
    def drij(self):
        """
        Element drij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 239
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__drij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            drij = self._arrays[array_handle]
        else:
            drij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__drij)
            self._arrays[array_handle] = drij
        return drij
    
    @drij.setter
    def drij(self, drij):
        self.drij[...] = drij
    
    @property
    def dzij(self):
        """
        Element dzij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 239
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzij = self._arrays[array_handle]
        else:
            dzij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzij)
            self._arrays[array_handle] = dzij
        return dzij
    
    @dzij.setter
    def dzij(self, dzij):
        self.dzij[...] = dzij
    
    @property
    def trij(self):
        """
        Element trij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 240
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__trij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            trij = self._arrays[array_handle]
        else:
            trij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__trij)
            self._arrays[array_handle] = trij
        return trij
    
    @trij.setter
    def trij(self, trij):
        self.trij[...] = trij
    
    @property
    def tzij(self):
        """
        Element tzij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 240
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tzij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tzij = self._arrays[array_handle]
        else:
            tzij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tzij)
            self._arrays[array_handle] = tzij
        return tzij
    
    @tzij.setter
    def tzij(self, tzij):
        self.tzij[...] = tzij
    
    @property
    def ivns(self):
        """
        Element ivns ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 241
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ivns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ivns = self._arrays[array_handle]
        else:
            ivns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ivns)
            self._arrays[array_handle] = ivns
        return ivns
    
    @ivns.setter
    def ivns(self, ivns):
        self.ivns[...] = ivns
    
    @property
    def ibns(self):
        """
        Element ibns ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 242
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ibns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ibns = self._arrays[array_handle]
        else:
            ibns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ibns)
            self._arrays[array_handle] = ibns
        return ibns
    
    @ibns.setter
    def ibns(self, ibns):
        self.ibns[...] = ibns
    
    @property
    def ivnc(self):
        """
        Element ivnc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 243
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ivnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ivnc = self._arrays[array_handle]
        else:
            ivnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ivnc)
            self._arrays[array_handle] = ivnc
        return ivnc
    
    @ivnc.setter
    def ivnc(self, ivnc):
        self.ivnc[...] = ivnc
    
    @property
    def ibnc(self):
        """
        Element ibnc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 244
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ibnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ibnc = self._arrays[array_handle]
        else:
            ibnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ibnc)
            self._arrays[array_handle] = ibnc
        return ibnc
    
    @ibnc.setter
    def ibnc(self, ibnc):
        self.ibnc[...] = ibnc
    
    @property
    def lrbc(self):
        """
        Element lrbc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 245
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lrbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lrbc = self._arrays[array_handle]
        else:
            lrbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lrbc)
            self._arrays[array_handle] = lrbc
        return lrbc
    
    @lrbc.setter
    def lrbc(self, lrbc):
        self.lrbc[...] = lrbc
    
    @property
    def lzbs(self):
        """
        Element lzbs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 245
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lzbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lzbs = self._arrays[array_handle]
        else:
            lzbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lzbs)
            self._arrays[array_handle] = lzbs
        return lzbs
    
    @lzbs.setter
    def lzbs(self, lzbs):
        self.lzbs[...] = lzbs
    
    @property
    def lrbs(self):
        """
        Element lrbs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 246
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lrbs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lrbs = self._arrays[array_handle]
        else:
            lrbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lrbs)
            self._arrays[array_handle] = lrbs
        return lrbs
    
    @lrbs.setter
    def lrbs(self, lrbs):
        self.lrbs[...] = lrbs
    
    @property
    def lzbc(self):
        """
        Element lzbc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 246
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lzbc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lzbc = self._arrays[array_handle]
        else:
            lzbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lzbc)
            self._arrays[array_handle] = lzbc
        return lzbc
    
    @lzbc.setter
    def lzbc(self, lzbc):
        self.lzbc[...] = lzbc
    
    @property
    def num_modes(self):
        """
        Element num_modes ftype=integer               pytype=int
        
        
        Defined at global.fpp line 248
        
        """
        return _spec.f90wrap_allglobal__get__num_modes()
    
    @num_modes.setter
    def num_modes(self, num_modes):
        _spec.f90wrap_allglobal__set__num_modes(num_modes)
    
    @property
    def mmrzrz(self):
        """
        Element mmrzrz ftype=integer pytype=int
        
        
        Defined at global.fpp line 249
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__mmrzrz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            mmrzrz = self._arrays[array_handle]
        else:
            mmrzrz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__mmrzrz)
            self._arrays[array_handle] = mmrzrz
        return mmrzrz
    
    @mmrzrz.setter
    def mmrzrz(self, mmrzrz):
        self.mmrzrz[...] = mmrzrz
    
    @property
    def nnrzrz(self):
        """
        Element nnrzrz ftype=integer pytype=int
        
        
        Defined at global.fpp line 249
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__nnrzrz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            nnrzrz = self._arrays[array_handle]
        else:
            nnrzrz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__nnrzrz)
            self._arrays[array_handle] = nnrzrz
        return nnrzrz
    
    @nnrzrz.setter
    def nnrzrz(self, nnrzrz):
        self.nnrzrz[...] = nnrzrz
    
    @property
    def allrzrz(self):
        """
        Element allrzrz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 250
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__allrzrz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            allrzrz = self._arrays[array_handle]
        else:
            allrzrz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__allrzrz)
            self._arrays[array_handle] = allrzrz
        return allrzrz
    
    @allrzrz.setter
    def allrzrz(self, allrzrz):
        self.allrzrz[...] = allrzrz
    
    @property
    def nt(self):
        """
        Element nt ftype=integer               pytype=int
        
        
        Defined at global.fpp line 256
        
        """
        return _spec.f90wrap_allglobal__get__nt()
    
    @nt.setter
    def nt(self, nt):
        _spec.f90wrap_allglobal__set__nt(nt)
    
    @property
    def nz(self):
        """
        Element nz ftype=integer               pytype=int
        
        
        Defined at global.fpp line 256
        
        """
        return _spec.f90wrap_allglobal__get__nz()
    
    @nz.setter
    def nz(self, nz):
        _spec.f90wrap_allglobal__set__nz(nz)
    
    @property
    def ntz(self):
        """
        Element ntz ftype=integer               pytype=int
        
        
        Defined at global.fpp line 256
        
        """
        return _spec.f90wrap_allglobal__get__ntz()
    
    @ntz.setter
    def ntz(self, ntz):
        _spec.f90wrap_allglobal__set__ntz(ntz)
    
    @property
    def hnt(self):
        """
        Element hnt ftype=integer               pytype=int
        
        
        Defined at global.fpp line 256
        
        """
        return _spec.f90wrap_allglobal__get__hnt()
    
    @hnt.setter
    def hnt(self, hnt):
        _spec.f90wrap_allglobal__set__hnt(hnt)
    
    @property
    def hnz(self):
        """
        Element hnz ftype=integer               pytype=int
        
        
        Defined at global.fpp line 256
        
        """
        return _spec.f90wrap_allglobal__get__hnz()
    
    @hnz.setter
    def hnz(self, hnz):
        _spec.f90wrap_allglobal__set__hnz(hnz)
    
    @property
    def sontz(self):
        """
        Element sontz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 257
        
        """
        return _spec.f90wrap_allglobal__get__sontz()
    
    @sontz.setter
    def sontz(self, sontz):
        _spec.f90wrap_allglobal__set__sontz(sontz)
    
    @property
    def rij(self):
        """
        Element rij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__rij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rij = self._arrays[array_handle]
        else:
            rij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__rij)
            self._arrays[array_handle] = rij
        return rij
    
    @rij.setter
    def rij(self, rij):
        self.rij[...] = rij
    
    @property
    def zij(self):
        """
        Element zij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__zij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zij = self._arrays[array_handle]
        else:
            zij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__zij)
            self._arrays[array_handle] = zij
        return zij
    
    @zij.setter
    def zij(self, zij):
        self.zij[...] = zij
    
    @property
    def xij(self):
        """
        Element xij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__xij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            xij = self._arrays[array_handle]
        else:
            xij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__xij)
            self._arrays[array_handle] = xij
        return xij
    
    @xij.setter
    def xij(self, xij):
        self.xij[...] = xij
    
    @property
    def yij(self):
        """
        Element yij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__yij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            yij = self._arrays[array_handle]
        else:
            yij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__yij)
            self._arrays[array_handle] = yij
        return yij
    
    @yij.setter
    def yij(self, yij):
        self.yij[...] = yij
    
    @property
    def sg(self):
        """
        Element sg ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__sg(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sg = self._arrays[array_handle]
        else:
            sg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__sg)
            self._arrays[array_handle] = sg
        return sg
    
    @sg.setter
    def sg(self, sg):
        self.sg[...] = sg
    
    @property
    def guvij(self):
        """
        Element guvij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__guvij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            guvij = self._arrays[array_handle]
        else:
            guvij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__guvij)
            self._arrays[array_handle] = guvij
        return guvij
    
    @guvij.setter
    def guvij(self, guvij):
        self.guvij[...] = guvij
    
    @property
    def gvuij(self):
        """
        Element gvuij ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 263
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gvuij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gvuij = self._arrays[array_handle]
        else:
            gvuij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gvuij)
            self._arrays[array_handle] = gvuij
        return gvuij
    
    @gvuij.setter
    def gvuij(self, gvuij):
        self.gvuij[...] = gvuij
    
    @property
    def guvijsave(self):
        """
        Element guvijsave ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 264
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__guvijsave(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            guvijsave = self._arrays[array_handle]
        else:
            guvijsave = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__guvijsave)
            self._arrays[array_handle] = guvijsave
        return guvijsave
    
    @guvijsave.setter
    def guvijsave(self, guvijsave):
        self.guvijsave[...] = guvijsave
    
    @property
    def ki(self):
        """
        Element ki ftype=integer pytype=int
        
        
        Defined at global.fpp line 265
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ki(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ki = self._arrays[array_handle]
        else:
            ki = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ki)
            self._arrays[array_handle] = ki
        return ki
    
    @ki.setter
    def ki(self, ki):
        self.ki[...] = ki
    
    @property
    def kijs(self):
        """
        Element kijs ftype=integer pytype=int
        
        
        Defined at global.fpp line 265
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__kijs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kijs = self._arrays[array_handle]
        else:
            kijs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__kijs)
            self._arrays[array_handle] = kijs
        return kijs
    
    @kijs.setter
    def kijs(self, kijs):
        self.kijs[...] = kijs
    
    @property
    def kija(self):
        """
        Element kija ftype=integer pytype=int
        
        
        Defined at global.fpp line 265
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__kija(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kija = self._arrays[array_handle]
        else:
            kija = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__kija)
            self._arrays[array_handle] = kija
        return kija
    
    @kija.setter
    def kija(self, kija):
        self.kija[...] = kija
    
    @property
    def iotakkii(self):
        """
        Element iotakkii ftype=integer pytype=int
        
        
        Defined at global.fpp line 266
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iotakkii(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iotakkii = self._arrays[array_handle]
        else:
            iotakkii = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iotakkii)
            self._arrays[array_handle] = iotakkii
        return iotakkii
    
    @iotakkii.setter
    def iotakkii(self, iotakkii):
        self.iotakkii[...] = iotakkii
    
    @property
    def iotaksub(self):
        """
        Element iotaksub ftype=integer pytype=int
        
        
        Defined at global.fpp line 266
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iotaksub(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iotaksub = self._arrays[array_handle]
        else:
            iotaksub = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iotaksub)
            self._arrays[array_handle] = iotaksub
        return iotaksub
    
    @iotaksub.setter
    def iotaksub(self, iotaksub):
        self.iotaksub[...] = iotaksub
    
    @property
    def iotakadd(self):
        """
        Element iotakadd ftype=integer pytype=int
        
        
        Defined at global.fpp line 266
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iotakadd(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iotakadd = self._arrays[array_handle]
        else:
            iotakadd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iotakadd)
            self._arrays[array_handle] = iotakadd
        return iotakadd
    
    @iotakadd.setter
    def iotakadd(self, iotakadd):
        self.iotakadd[...] = iotakadd
    
    @property
    def iotaksgn(self):
        """
        Element iotaksgn ftype=integer pytype=int
        
        
        Defined at global.fpp line 266
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iotaksgn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iotaksgn = self._arrays[array_handle]
        else:
            iotaksgn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iotaksgn)
            self._arrays[array_handle] = iotaksgn
        return iotaksgn
    
    @iotaksgn.setter
    def iotaksgn(self, iotaksgn):
        self.iotaksgn[...] = iotaksgn
    
    @property
    def efmn(self):
        """
        Element efmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 267
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__efmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            efmn = self._arrays[array_handle]
        else:
            efmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__efmn)
            self._arrays[array_handle] = efmn
        return efmn
    
    @efmn.setter
    def efmn(self, efmn):
        self.efmn[...] = efmn
    
    @property
    def ofmn(self):
        """
        Element ofmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 267
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ofmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ofmn = self._arrays[array_handle]
        else:
            ofmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ofmn)
            self._arrays[array_handle] = ofmn
        return ofmn
    
    @ofmn.setter
    def ofmn(self, ofmn):
        self.ofmn[...] = ofmn
    
    @property
    def cfmn(self):
        """
        Element cfmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 267
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__cfmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cfmn = self._arrays[array_handle]
        else:
            cfmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__cfmn)
            self._arrays[array_handle] = cfmn
        return cfmn
    
    @cfmn.setter
    def cfmn(self, cfmn):
        self.cfmn[...] = cfmn
    
    @property
    def sfmn(self):
        """
        Element sfmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 267
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__sfmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sfmn = self._arrays[array_handle]
        else:
            sfmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__sfmn)
            self._arrays[array_handle] = sfmn
        return sfmn
    
    @sfmn.setter
    def sfmn(self, sfmn):
        self.sfmn[...] = sfmn
    
    @property
    def evmn(self):
        """
        Element evmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 268
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__evmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            evmn = self._arrays[array_handle]
        else:
            evmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__evmn)
            self._arrays[array_handle] = evmn
        return evmn
    
    @evmn.setter
    def evmn(self, evmn):
        self.evmn[...] = evmn
    
    @property
    def odmn(self):
        """
        Element odmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 268
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__odmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            odmn = self._arrays[array_handle]
        else:
            odmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__odmn)
            self._arrays[array_handle] = odmn
        return odmn
    
    @odmn.setter
    def odmn(self, odmn):
        self.odmn[...] = odmn
    
    @property
    def comn(self):
        """
        Element comn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 268
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__comn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            comn = self._arrays[array_handle]
        else:
            comn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__comn)
            self._arrays[array_handle] = comn
        return comn
    
    @comn.setter
    def comn(self, comn):
        self.comn[...] = comn
    
    @property
    def simn(self):
        """
        Element simn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 268
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__simn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            simn = self._arrays[array_handle]
        else:
            simn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__simn)
            self._arrays[array_handle] = simn
        return simn
    
    @simn.setter
    def simn(self, simn):
        self.simn[...] = simn
    
    @property
    def ijreal(self):
        """
        Element ijreal ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 269
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ijreal(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ijreal = self._arrays[array_handle]
        else:
            ijreal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ijreal)
            self._arrays[array_handle] = ijreal
        return ijreal
    
    @ijreal.setter
    def ijreal(self, ijreal):
        self.ijreal[...] = ijreal
    
    @property
    def ijimag(self):
        """
        Element ijimag ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 269
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ijimag(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ijimag = self._arrays[array_handle]
        else:
            ijimag = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ijimag)
            self._arrays[array_handle] = ijimag
        return ijimag
    
    @ijimag.setter
    def ijimag(self, ijimag):
        self.ijimag[...] = ijimag
    
    @property
    def jireal(self):
        """
        Element jireal ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 269
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__jireal(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jireal = self._arrays[array_handle]
        else:
            jireal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__jireal)
            self._arrays[array_handle] = jireal
        return jireal
    
    @jireal.setter
    def jireal(self, jireal):
        self.jireal[...] = jireal
    
    @property
    def jiimag(self):
        """
        Element jiimag ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 269
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__jiimag(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jiimag = self._arrays[array_handle]
        else:
            jiimag = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__jiimag)
            self._arrays[array_handle] = jiimag
        return jiimag
    
    @jiimag.setter
    def jiimag(self, jiimag):
        self.jiimag[...] = jiimag
    
    @property
    def jkreal(self):
        """
        Element jkreal ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 270
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__jkreal(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jkreal = self._arrays[array_handle]
        else:
            jkreal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__jkreal)
            self._arrays[array_handle] = jkreal
        return jkreal
    
    @jkreal.setter
    def jkreal(self, jkreal):
        self.jkreal[...] = jkreal
    
    @property
    def jkimag(self):
        """
        Element jkimag ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 270
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__jkimag(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jkimag = self._arrays[array_handle]
        else:
            jkimag = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__jkimag)
            self._arrays[array_handle] = jkimag
        return jkimag
    
    @jkimag.setter
    def jkimag(self, jkimag):
        self.jkimag[...] = jkimag
    
    @property
    def kjreal(self):
        """
        Element kjreal ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 270
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__kjreal(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kjreal = self._arrays[array_handle]
        else:
            kjreal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__kjreal)
            self._arrays[array_handle] = kjreal
        return kjreal
    
    @kjreal.setter
    def kjreal(self, kjreal):
        self.kjreal[...] = kjreal
    
    @property
    def kjimag(self):
        """
        Element kjimag ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 270
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__kjimag(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kjimag = self._arrays[array_handle]
        else:
            kjimag = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__kjimag)
            self._arrays[array_handle] = kjimag
        return kjimag
    
    @kjimag.setter
    def kjimag(self, kjimag):
        self.kjimag[...] = kjimag
    
    @property
    def bsupumn(self):
        """
        Element bsupumn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 271
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bsupumn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bsupumn = self._arrays[array_handle]
        else:
            bsupumn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bsupumn)
            self._arrays[array_handle] = bsupumn
        return bsupumn
    
    @bsupumn.setter
    def bsupumn(self, bsupumn):
        self.bsupumn[...] = bsupumn
    
    @property
    def bsupvmn(self):
        """
        Element bsupvmn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 271
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bsupvmn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bsupvmn = self._arrays[array_handle]
        else:
            bsupvmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bsupvmn)
            self._arrays[array_handle] = bsupvmn
        return bsupvmn
    
    @bsupvmn.setter
    def bsupvmn(self, bsupvmn):
        self.bsupvmn[...] = bsupvmn
    
    @property
    def goomne(self):
        """
        Element goomne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 273
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__goomne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            goomne = self._arrays[array_handle]
        else:
            goomne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__goomne)
            self._arrays[array_handle] = goomne
        return goomne
    
    @goomne.setter
    def goomne(self, goomne):
        self.goomne[...] = goomne
    
    @property
    def goomno(self):
        """
        Element goomno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 273
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__goomno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            goomno = self._arrays[array_handle]
        else:
            goomno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__goomno)
            self._arrays[array_handle] = goomno
        return goomno
    
    @goomno.setter
    def goomno(self, goomno):
        self.goomno[...] = goomno
    
    @property
    def gssmne(self):
        """
        Element gssmne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 274
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gssmne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gssmne = self._arrays[array_handle]
        else:
            gssmne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gssmne)
            self._arrays[array_handle] = gssmne
        return gssmne
    
    @gssmne.setter
    def gssmne(self, gssmne):
        self.gssmne[...] = gssmne
    
    @property
    def gssmno(self):
        """
        Element gssmno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 274
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gssmno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gssmno = self._arrays[array_handle]
        else:
            gssmno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gssmno)
            self._arrays[array_handle] = gssmno
        return gssmno
    
    @gssmno.setter
    def gssmno(self, gssmno):
        self.gssmno[...] = gssmno
    
    @property
    def gstmne(self):
        """
        Element gstmne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 275
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gstmne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gstmne = self._arrays[array_handle]
        else:
            gstmne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gstmne)
            self._arrays[array_handle] = gstmne
        return gstmne
    
    @gstmne.setter
    def gstmne(self, gstmne):
        self.gstmne[...] = gstmne
    
    @property
    def gstmno(self):
        """
        Element gstmno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 275
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gstmno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gstmno = self._arrays[array_handle]
        else:
            gstmno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gstmno)
            self._arrays[array_handle] = gstmno
        return gstmno
    
    @gstmno.setter
    def gstmno(self, gstmno):
        self.gstmno[...] = gstmno
    
    @property
    def gszmne(self):
        """
        Element gszmne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 276
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gszmne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gszmne = self._arrays[array_handle]
        else:
            gszmne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gszmne)
            self._arrays[array_handle] = gszmne
        return gszmne
    
    @gszmne.setter
    def gszmne(self, gszmne):
        self.gszmne[...] = gszmne
    
    @property
    def gszmno(self):
        """
        Element gszmno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 276
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gszmno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gszmno = self._arrays[array_handle]
        else:
            gszmno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gszmno)
            self._arrays[array_handle] = gszmno
        return gszmno
    
    @gszmno.setter
    def gszmno(self, gszmno):
        self.gszmno[...] = gszmno
    
    @property
    def gttmne(self):
        """
        Element gttmne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 277
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gttmne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gttmne = self._arrays[array_handle]
        else:
            gttmne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gttmne)
            self._arrays[array_handle] = gttmne
        return gttmne
    
    @gttmne.setter
    def gttmne(self, gttmne):
        self.gttmne[...] = gttmne
    
    @property
    def gttmno(self):
        """
        Element gttmno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 277
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gttmno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gttmno = self._arrays[array_handle]
        else:
            gttmno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gttmno)
            self._arrays[array_handle] = gttmno
        return gttmno
    
    @gttmno.setter
    def gttmno(self, gttmno):
        self.gttmno[...] = gttmno
    
    @property
    def gtzmne(self):
        """
        Element gtzmne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 278
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gtzmne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gtzmne = self._arrays[array_handle]
        else:
            gtzmne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gtzmne)
            self._arrays[array_handle] = gtzmne
        return gtzmne
    
    @gtzmne.setter
    def gtzmne(self, gtzmne):
        self.gtzmne[...] = gtzmne
    
    @property
    def gtzmno(self):
        """
        Element gtzmno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 278
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gtzmno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gtzmno = self._arrays[array_handle]
        else:
            gtzmno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gtzmno)
            self._arrays[array_handle] = gtzmno
        return gtzmno
    
    @gtzmno.setter
    def gtzmno(self, gtzmno):
        self.gtzmno[...] = gtzmno
    
    @property
    def gzzmne(self):
        """
        Element gzzmne ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 279
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gzzmne(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gzzmne = self._arrays[array_handle]
        else:
            gzzmne = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gzzmne)
            self._arrays[array_handle] = gzzmne
        return gzzmne
    
    @gzzmne.setter
    def gzzmne(self, gzzmne):
        self.gzzmne[...] = gzzmne
    
    @property
    def gzzmno(self):
        """
        Element gzzmno ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 279
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gzzmno(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gzzmno = self._arrays[array_handle]
        else:
            gzzmno = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gzzmno)
            self._arrays[array_handle] = gzzmno
        return gzzmno
    
    @gzzmno.setter
    def gzzmno(self, gzzmno):
        self.gzzmno[...] = gzzmno
    
    @property
    def dtoocc(self):
        """
        Element dtoocc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 291
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dtoocc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dtoocc = self._arrays[array_handle]
        else:
            dtoocc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dtoocc)
            self._arrays[array_handle] = dtoocc
        return dtoocc
    
    @dtoocc.setter
    def dtoocc(self, dtoocc):
        self.dtoocc[...] = dtoocc
    
    @property
    def dtoocs(self):
        """
        Element dtoocs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 291
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dtoocs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dtoocs = self._arrays[array_handle]
        else:
            dtoocs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dtoocs)
            self._arrays[array_handle] = dtoocs
        return dtoocs
    
    @dtoocs.setter
    def dtoocs(self, dtoocs):
        self.dtoocs[...] = dtoocs
    
    @property
    def dtoosc(self):
        """
        Element dtoosc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 291
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dtoosc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dtoosc = self._arrays[array_handle]
        else:
            dtoosc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dtoosc)
            self._arrays[array_handle] = dtoosc
        return dtoosc
    
    @dtoosc.setter
    def dtoosc(self, dtoosc):
        self.dtoosc[...] = dtoosc
    
    @property
    def dtooss(self):
        """
        Element dtooss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 291
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dtooss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dtooss = self._arrays[array_handle]
        else:
            dtooss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dtooss)
            self._arrays[array_handle] = dtooss
        return dtooss
    
    @dtooss.setter
    def dtooss(self, dtooss):
        self.dtooss[...] = dtooss
    
    @property
    def ttsscc(self):
        """
        Element ttsscc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 292
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ttsscc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ttsscc = self._arrays[array_handle]
        else:
            ttsscc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ttsscc)
            self._arrays[array_handle] = ttsscc
        return ttsscc
    
    @ttsscc.setter
    def ttsscc(self, ttsscc):
        self.ttsscc[...] = ttsscc
    
    @property
    def ttsscs(self):
        """
        Element ttsscs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 292
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ttsscs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ttsscs = self._arrays[array_handle]
        else:
            ttsscs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ttsscs)
            self._arrays[array_handle] = ttsscs
        return ttsscs
    
    @ttsscs.setter
    def ttsscs(self, ttsscs):
        self.ttsscs[...] = ttsscs
    
    @property
    def ttsssc(self):
        """
        Element ttsssc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 292
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ttsssc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ttsssc = self._arrays[array_handle]
        else:
            ttsssc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ttsssc)
            self._arrays[array_handle] = ttsssc
        return ttsssc
    
    @ttsssc.setter
    def ttsssc(self, ttsssc):
        self.ttsssc[...] = ttsssc
    
    @property
    def ttssss(self):
        """
        Element ttssss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 292
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ttssss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ttssss = self._arrays[array_handle]
        else:
            ttssss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ttssss)
            self._arrays[array_handle] = ttssss
        return ttssss
    
    @ttssss.setter
    def ttssss(self, ttssss):
        self.ttssss[...] = ttssss
    
    @property
    def tdstcc(self):
        """
        Element tdstcc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 293
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdstcc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdstcc = self._arrays[array_handle]
        else:
            tdstcc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdstcc)
            self._arrays[array_handle] = tdstcc
        return tdstcc
    
    @tdstcc.setter
    def tdstcc(self, tdstcc):
        self.tdstcc[...] = tdstcc
    
    @property
    def tdstcs(self):
        """
        Element tdstcs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 293
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdstcs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdstcs = self._arrays[array_handle]
        else:
            tdstcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdstcs)
            self._arrays[array_handle] = tdstcs
        return tdstcs
    
    @tdstcs.setter
    def tdstcs(self, tdstcs):
        self.tdstcs[...] = tdstcs
    
    @property
    def tdstsc(self):
        """
        Element tdstsc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 293
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdstsc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdstsc = self._arrays[array_handle]
        else:
            tdstsc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdstsc)
            self._arrays[array_handle] = tdstsc
        return tdstsc
    
    @tdstsc.setter
    def tdstsc(self, tdstsc):
        self.tdstsc[...] = tdstsc
    
    @property
    def tdstss(self):
        """
        Element tdstss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 293
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdstss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdstss = self._arrays[array_handle]
        else:
            tdstss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdstss)
            self._arrays[array_handle] = tdstss
        return tdstss
    
    @tdstss.setter
    def tdstss(self, tdstss):
        self.tdstss[...] = tdstss
    
    @property
    def tdszcc(self):
        """
        Element tdszcc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 294
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdszcc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdszcc = self._arrays[array_handle]
        else:
            tdszcc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdszcc)
            self._arrays[array_handle] = tdszcc
        return tdszcc
    
    @tdszcc.setter
    def tdszcc(self, tdszcc):
        self.tdszcc[...] = tdszcc
    
    @property
    def tdszcs(self):
        """
        Element tdszcs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 294
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdszcs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdszcs = self._arrays[array_handle]
        else:
            tdszcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdszcs)
            self._arrays[array_handle] = tdszcs
        return tdszcs
    
    @tdszcs.setter
    def tdszcs(self, tdszcs):
        self.tdszcs[...] = tdszcs
    
    @property
    def tdszsc(self):
        """
        Element tdszsc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 294
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdszsc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdszsc = self._arrays[array_handle]
        else:
            tdszsc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdszsc)
            self._arrays[array_handle] = tdszsc
        return tdszsc
    
    @tdszsc.setter
    def tdszsc(self, tdszsc):
        self.tdszsc[...] = tdszsc
    
    @property
    def tdszss(self):
        """
        Element tdszss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 294
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tdszss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tdszss = self._arrays[array_handle]
        else:
            tdszss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tdszss)
            self._arrays[array_handle] = tdszss
        return tdszss
    
    @tdszss.setter
    def tdszss(self, tdszss):
        self.tdszss[...] = tdszss
    
    @property
    def ddttcc(self):
        """
        Element ddttcc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 295
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddttcc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddttcc = self._arrays[array_handle]
        else:
            ddttcc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddttcc)
            self._arrays[array_handle] = ddttcc
        return ddttcc
    
    @ddttcc.setter
    def ddttcc(self, ddttcc):
        self.ddttcc[...] = ddttcc
    
    @property
    def ddttcs(self):
        """
        Element ddttcs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 295
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddttcs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddttcs = self._arrays[array_handle]
        else:
            ddttcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddttcs)
            self._arrays[array_handle] = ddttcs
        return ddttcs
    
    @ddttcs.setter
    def ddttcs(self, ddttcs):
        self.ddttcs[...] = ddttcs
    
    @property
    def ddttsc(self):
        """
        Element ddttsc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 295
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddttsc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddttsc = self._arrays[array_handle]
        else:
            ddttsc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddttsc)
            self._arrays[array_handle] = ddttsc
        return ddttsc
    
    @ddttsc.setter
    def ddttsc(self, ddttsc):
        self.ddttsc[...] = ddttsc
    
    @property
    def ddttss(self):
        """
        Element ddttss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 295
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddttss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddttss = self._arrays[array_handle]
        else:
            ddttss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddttss)
            self._arrays[array_handle] = ddttss
        return ddttss
    
    @ddttss.setter
    def ddttss(self, ddttss):
        self.ddttss[...] = ddttss
    
    @property
    def ddtzcc(self):
        """
        Element ddtzcc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 296
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddtzcc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddtzcc = self._arrays[array_handle]
        else:
            ddtzcc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddtzcc)
            self._arrays[array_handle] = ddtzcc
        return ddtzcc
    
    @ddtzcc.setter
    def ddtzcc(self, ddtzcc):
        self.ddtzcc[...] = ddtzcc
    
    @property
    def ddtzcs(self):
        """
        Element ddtzcs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 296
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddtzcs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddtzcs = self._arrays[array_handle]
        else:
            ddtzcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddtzcs)
            self._arrays[array_handle] = ddtzcs
        return ddtzcs
    
    @ddtzcs.setter
    def ddtzcs(self, ddtzcs):
        self.ddtzcs[...] = ddtzcs
    
    @property
    def ddtzsc(self):
        """
        Element ddtzsc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 296
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddtzsc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddtzsc = self._arrays[array_handle]
        else:
            ddtzsc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddtzsc)
            self._arrays[array_handle] = ddtzsc
        return ddtzsc
    
    @ddtzsc.setter
    def ddtzsc(self, ddtzsc):
        self.ddtzsc[...] = ddtzsc
    
    @property
    def ddtzss(self):
        """
        Element ddtzss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 296
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddtzss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddtzss = self._arrays[array_handle]
        else:
            ddtzss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddtzss)
            self._arrays[array_handle] = ddtzss
        return ddtzss
    
    @ddtzss.setter
    def ddtzss(self, ddtzss):
        self.ddtzss[...] = ddtzss
    
    @property
    def ddzzcc(self):
        """
        Element ddzzcc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 297
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddzzcc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddzzcc = self._arrays[array_handle]
        else:
            ddzzcc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddzzcc)
            self._arrays[array_handle] = ddzzcc
        return ddzzcc
    
    @ddzzcc.setter
    def ddzzcc(self, ddzzcc):
        self.ddzzcc[...] = ddzzcc
    
    @property
    def ddzzcs(self):
        """
        Element ddzzcs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 297
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddzzcs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddzzcs = self._arrays[array_handle]
        else:
            ddzzcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddzzcs)
            self._arrays[array_handle] = ddzzcs
        return ddzzcs
    
    @ddzzcs.setter
    def ddzzcs(self, ddzzcs):
        self.ddzzcs[...] = ddzzcs
    
    @property
    def ddzzsc(self):
        """
        Element ddzzsc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 297
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddzzsc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddzzsc = self._arrays[array_handle]
        else:
            ddzzsc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddzzsc)
            self._arrays[array_handle] = ddzzsc
        return ddzzsc
    
    @ddzzsc.setter
    def ddzzsc(self, ddzzsc):
        self.ddzzsc[...] = ddzzsc
    
    @property
    def ddzzss(self):
        """
        Element ddzzss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 297
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ddzzss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ddzzss = self._arrays[array_handle]
        else:
            ddzzss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ddzzss)
            self._arrays[array_handle] = ddzzss
        return ddzzss
    
    @ddzzss.setter
    def ddzzss(self, ddzzss):
        self.ddzzss[...] = ddzzss
    
    @property
    def tsc(self):
        """
        Element tsc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 299
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tsc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tsc = self._arrays[array_handle]
        else:
            tsc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tsc)
            self._arrays[array_handle] = tsc
        return tsc
    
    @tsc.setter
    def tsc(self, tsc):
        self.tsc[...] = tsc
    
    @property
    def tss(self):
        """
        Element tss ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 299
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tss(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tss = self._arrays[array_handle]
        else:
            tss = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tss)
            self._arrays[array_handle] = tss
        return tss
    
    @tss.setter
    def tss(self, tss):
        self.tss[...] = tss
    
    @property
    def dtc(self):
        """
        Element dtc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 299
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dtc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dtc = self._arrays[array_handle]
        else:
            dtc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dtc)
            self._arrays[array_handle] = dtc
        return dtc
    
    @dtc.setter
    def dtc(self, dtc):
        self.dtc[...] = dtc
    
    @property
    def dts(self):
        """
        Element dts ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 299
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dts(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dts = self._arrays[array_handle]
        else:
            dts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dts)
            self._arrays[array_handle] = dts
        return dts
    
    @dts.setter
    def dts(self, dts):
        self.dts[...] = dts
    
    @property
    def dzc(self):
        """
        Element dzc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 299
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzc = self._arrays[array_handle]
        else:
            dzc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzc)
            self._arrays[array_handle] = dzc
        return dzc
    
    @dzc.setter
    def dzc(self, dzc):
        self.dzc[...] = dzc
    
    @property
    def dzs(self):
        """
        Element dzs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 299
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzs = self._arrays[array_handle]
        else:
            dzs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzs)
            self._arrays[array_handle] = dzs
        return dzs
    
    @dzs.setter
    def dzs(self, dzs):
        self.dzs[...] = dzs
    
    @property
    def ttc(self):
        """
        Element ttc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 300
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ttc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ttc = self._arrays[array_handle]
        else:
            ttc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ttc)
            self._arrays[array_handle] = ttc
        return ttc
    
    @ttc.setter
    def ttc(self, ttc):
        self.ttc[...] = ttc
    
    @property
    def tzc(self):
        """
        Element tzc ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 300
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tzc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tzc = self._arrays[array_handle]
        else:
            tzc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tzc)
            self._arrays[array_handle] = tzc
        return tzc
    
    @tzc.setter
    def tzc(self, tzc):
        self.tzc[...] = tzc
    
    @property
    def tts(self):
        """
        Element tts ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 300
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tts(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tts = self._arrays[array_handle]
        else:
            tts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tts)
            self._arrays[array_handle] = tts
        return tts
    
    @tts.setter
    def tts(self, tts):
        self.tts[...] = tts
    
    @property
    def tzs(self):
        """
        Element tzs ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 300
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tzs(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tzs = self._arrays[array_handle]
        else:
            tzs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tzs)
            self._arrays[array_handle] = tzs
        return tzs
    
    @tzs.setter
    def tzs(self, tzs):
        self.tzs[...] = tzs
    
    @property
    def dtflux(self):
        """
        Element dtflux ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 302
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dtflux(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dtflux = self._arrays[array_handle]
        else:
            dtflux = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dtflux)
            self._arrays[array_handle] = dtflux
        return dtflux
    
    @dtflux.setter
    def dtflux(self, dtflux):
        self.dtflux[...] = dtflux
    
    @property
    def dpflux(self):
        """
        Element dpflux ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 302
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dpflux(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dpflux = self._arrays[array_handle]
        else:
            dpflux = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dpflux)
            self._arrays[array_handle] = dpflux
        return dpflux
    
    @dpflux.setter
    def dpflux(self, dpflux):
        self.dpflux[...] = dpflux
    
    @property
    def sweight(self):
        """
        Element sweight ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 304
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__sweight(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sweight = self._arrays[array_handle]
        else:
            sweight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__sweight)
            self._arrays[array_handle] = sweight
        return sweight
    
    @sweight.setter
    def sweight(self, sweight):
        self.sweight[...] = sweight
    
    @property
    def nadof(self):
        """
        Element nadof ftype=integer pytype=int
        
        
        Defined at global.fpp line 310
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__nadof(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            nadof = self._arrays[array_handle]
        else:
            nadof = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__nadof)
            self._arrays[array_handle] = nadof
        return nadof
    
    @nadof.setter
    def nadof(self, nadof):
        self.nadof[...] = nadof
    
    @property
    def nfielddof(self):
        """
        Element nfielddof ftype=integer pytype=int
        
        
        Defined at global.fpp line 311
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__nfielddof(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            nfielddof = self._arrays[array_handle]
        else:
            nfielddof = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__nfielddof)
            self._arrays[array_handle] = nfielddof
        return nfielddof
    
    @nfielddof.setter
    def nfielddof(self, nfielddof):
        self.nfielddof[...] = nfielddof
    
    @property
    def lma(self):
        """
        Element lma ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lma(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lma = self._arrays[array_handle]
        else:
            lma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lma)
            self._arrays[array_handle] = lma
        return lma
    
    @lma.setter
    def lma(self, lma):
        self.lma[...] = lma
    
    @property
    def lmb(self):
        """
        Element lmb ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmb(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmb = self._arrays[array_handle]
        else:
            lmb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmb)
            self._arrays[array_handle] = lmb
        return lmb
    
    @lmb.setter
    def lmb(self, lmb):
        self.lmb[...] = lmb
    
    @property
    def lmc(self):
        """
        Element lmc ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmc = self._arrays[array_handle]
        else:
            lmc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmc)
            self._arrays[array_handle] = lmc
        return lmc
    
    @lmc.setter
    def lmc(self, lmc):
        self.lmc[...] = lmc
    
    @property
    def lmd(self):
        """
        Element lmd ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmd(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmd = self._arrays[array_handle]
        else:
            lmd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmd)
            self._arrays[array_handle] = lmd
        return lmd
    
    @lmd.setter
    def lmd(self, lmd):
        self.lmd[...] = lmd
    
    @property
    def lme(self):
        """
        Element lme ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lme(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lme = self._arrays[array_handle]
        else:
            lme = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lme)
            self._arrays[array_handle] = lme
        return lme
    
    @lme.setter
    def lme(self, lme):
        self.lme[...] = lme
    
    @property
    def lmf(self):
        """
        Element lmf ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmf(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmf = self._arrays[array_handle]
        else:
            lmf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmf)
            self._arrays[array_handle] = lmf
        return lmf
    
    @lmf.setter
    def lmf(self, lmf):
        self.lmf[...] = lmf
    
    @property
    def lmg(self):
        """
        Element lmg ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmg(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmg = self._arrays[array_handle]
        else:
            lmg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmg)
            self._arrays[array_handle] = lmg
        return lmg
    
    @lmg.setter
    def lmg(self, lmg):
        self.lmg[...] = lmg
    
    @property
    def lmh(self):
        """
        Element lmh ftype=integer       pytype=int
        
        
        Defined at global.fpp line 328
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmh(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmh = self._arrays[array_handle]
        else:
            lmh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmh)
            self._arrays[array_handle] = lmh
        return lmh
    
    @lmh.setter
    def lmh(self, lmh):
        self.lmh[...] = lmh
    
    @property
    def lmavalue(self):
        """
        Element lmavalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 329
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmavalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmavalue = self._arrays[array_handle]
        else:
            lmavalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmavalue)
            self._arrays[array_handle] = lmavalue
        return lmavalue
    
    @lmavalue.setter
    def lmavalue(self, lmavalue):
        self.lmavalue[...] = lmavalue
    
    @property
    def lmbvalue(self):
        """
        Element lmbvalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 329
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmbvalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmbvalue = self._arrays[array_handle]
        else:
            lmbvalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmbvalue)
            self._arrays[array_handle] = lmbvalue
        return lmbvalue
    
    @lmbvalue.setter
    def lmbvalue(self, lmbvalue):
        self.lmbvalue[...] = lmbvalue
    
    @property
    def lmcvalue(self):
        """
        Element lmcvalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 329
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmcvalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmcvalue = self._arrays[array_handle]
        else:
            lmcvalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmcvalue)
            self._arrays[array_handle] = lmcvalue
        return lmcvalue
    
    @lmcvalue.setter
    def lmcvalue(self, lmcvalue):
        self.lmcvalue[...] = lmcvalue
    
    @property
    def lmdvalue(self):
        """
        Element lmdvalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 329
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmdvalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmdvalue = self._arrays[array_handle]
        else:
            lmdvalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmdvalue)
            self._arrays[array_handle] = lmdvalue
        return lmdvalue
    
    @lmdvalue.setter
    def lmdvalue(self, lmdvalue):
        self.lmdvalue[...] = lmdvalue
    
    @property
    def lmevalue(self):
        """
        Element lmevalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 329
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmevalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmevalue = self._arrays[array_handle]
        else:
            lmevalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmevalue)
            self._arrays[array_handle] = lmevalue
        return lmevalue
    
    @lmevalue.setter
    def lmevalue(self, lmevalue):
        self.lmevalue[...] = lmevalue
    
    @property
    def lmfvalue(self):
        """
        Element lmfvalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 329
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmfvalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmfvalue = self._arrays[array_handle]
        else:
            lmfvalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmfvalue)
            self._arrays[array_handle] = lmfvalue
        return lmfvalue
    
    @lmfvalue.setter
    def lmfvalue(self, lmfvalue):
        self.lmfvalue[...] = lmfvalue
    
    @property
    def lmgvalue(self):
        """
        Element lmgvalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 330
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmgvalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmgvalue = self._arrays[array_handle]
        else:
            lmgvalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmgvalue)
            self._arrays[array_handle] = lmgvalue
        return lmgvalue
    
    @lmgvalue.setter
    def lmgvalue(self, lmgvalue):
        self.lmgvalue[...] = lmgvalue
    
    @property
    def lmhvalue(self):
        """
        Element lmhvalue ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 330
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lmhvalue(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lmhvalue = self._arrays[array_handle]
        else:
            lmhvalue = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lmhvalue)
            self._arrays[array_handle] = lmhvalue
        return lmhvalue
    
    @lmhvalue.setter
    def lmhvalue(self, lmhvalue):
        self.lmhvalue[...] = lmhvalue
    
    @property
    def fso(self):
        """
        Element fso ftype=integer       pytype=int
        
        
        Defined at global.fpp line 333
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__fso(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fso = self._arrays[array_handle]
        else:
            fso = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__fso)
            self._arrays[array_handle] = fso
        return fso
    
    @fso.setter
    def fso(self, fso):
        self.fso[...] = fso
    
    @property
    def fse(self):
        """
        Element fse ftype=integer       pytype=int
        
        
        Defined at global.fpp line 333
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__fse(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fse = self._arrays[array_handle]
        else:
            fse = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__fse)
            self._arrays[array_handle] = fse
        return fse
    
    @fse.setter
    def fse(self, fse):
        self.fse[...] = fse
    
    @property
    def lcoordinatesingularity(self):
        """
        Element lcoordinatesingularity ftype=logical pytype=bool
        
        
        Defined at global.fpp line 334
        
        """
        return _spec.f90wrap_allglobal__get__lcoordinatesingularity()
    
    @lcoordinatesingularity.setter
    def lcoordinatesingularity(self, lcoordinatesingularity):
        _spec.f90wrap_allglobal__set__lcoordinatesingularity(lcoordinatesingularity)
    
    @property
    def lplasmaregion(self):
        """
        Element lplasmaregion ftype=logical pytype=bool
        
        
        Defined at global.fpp line 334
        
        """
        return _spec.f90wrap_allglobal__get__lplasmaregion()
    
    @lplasmaregion.setter
    def lplasmaregion(self, lplasmaregion):
        _spec.f90wrap_allglobal__set__lplasmaregion(lplasmaregion)
    
    @property
    def lvacuumregion(self):
        """
        Element lvacuumregion ftype=logical pytype=bool
        
        
        Defined at global.fpp line 334
        
        """
        return _spec.f90wrap_allglobal__get__lvacuumregion()
    
    @lvacuumregion.setter
    def lvacuumregion(self, lvacuumregion):
        _spec.f90wrap_allglobal__set__lvacuumregion(lvacuumregion)
    
    @property
    def lsavedguvij(self):
        """
        Element lsavedguvij ftype=logical pytype=bool
        
        
        Defined at global.fpp line 335
        
        """
        return _spec.f90wrap_allglobal__get__lsavedguvij()
    
    @lsavedguvij.setter
    def lsavedguvij(self, lsavedguvij):
        _spec.f90wrap_allglobal__set__lsavedguvij(lsavedguvij)
    
    @property
    def localconstraint(self):
        """
        Element localconstraint ftype=logical pytype=bool
        
        
        Defined at global.fpp line 336
        
        """
        return _spec.f90wrap_allglobal__get__localconstraint()
    
    @localconstraint.setter
    def localconstraint(self, localconstraint):
        _spec.f90wrap_allglobal__set__localconstraint(localconstraint)
    
    @property
    def dma(self):
        """
        Element dma ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 347
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dma(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dma = self._arrays[array_handle]
        else:
            dma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dma)
            self._arrays[array_handle] = dma
        return dma
    
    @dma.setter
    def dma(self, dma):
        self.dma[...] = dma
    
    @property
    def dmb(self):
        """
        Element dmb ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 347
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dmb(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dmb = self._arrays[array_handle]
        else:
            dmb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dmb)
            self._arrays[array_handle] = dmb
        return dmb
    
    @dmb.setter
    def dmb(self, dmb):
        self.dmb[...] = dmb
    
    @property
    def dmd(self):
        """
        Element dmd ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 348
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dmd(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dmd = self._arrays[array_handle]
        else:
            dmd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dmd)
            self._arrays[array_handle] = dmd
        return dmd
    
    @dmd.setter
    def dmd(self, dmd):
        self.dmd[...] = dmd
    
    @property
    def dmas(self):
        """
        Element dmas ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 349
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dmas(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dmas = self._arrays[array_handle]
        else:
            dmas = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dmas)
            self._arrays[array_handle] = dmas
        return dmas
    
    @dmas.setter
    def dmas(self, dmas):
        self.dmas[...] = dmas
    
    @property
    def dmds(self):
        """
        Element dmds ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 349
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dmds(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dmds = self._arrays[array_handle]
        else:
            dmds = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dmds)
            self._arrays[array_handle] = dmds
        return dmds
    
    @dmds.setter
    def dmds(self, dmds):
        self.dmds[...] = dmds
    
    @property
    def idmas(self):
        """
        Element idmas ftype=integer pytype=int
        
        
        Defined at global.fpp line 350
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__idmas(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            idmas = self._arrays[array_handle]
        else:
            idmas = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__idmas)
            self._arrays[array_handle] = idmas
        return idmas
    
    @idmas.setter
    def idmas(self, idmas):
        self.idmas[...] = idmas
    
    @property
    def jdmas(self):
        """
        Element jdmas ftype=integer pytype=int
        
        
        Defined at global.fpp line 350
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__jdmas(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jdmas = self._arrays[array_handle]
        else:
            jdmas = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__jdmas)
            self._arrays[array_handle] = jdmas
        return jdmas
    
    @jdmas.setter
    def jdmas(self, jdmas):
        self.jdmas[...] = jdmas
    
    @property
    def ndmasmax(self):
        """
        Element ndmasmax ftype=integer pytype=int
        
        
        Defined at global.fpp line 351
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ndmasmax(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ndmasmax = self._arrays[array_handle]
        else:
            ndmasmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ndmasmax)
            self._arrays[array_handle] = ndmasmax
        return ndmasmax
    
    @ndmasmax.setter
    def ndmasmax(self, ndmasmax):
        self.ndmasmax[...] = ndmasmax
    
    @property
    def ndmas(self):
        """
        Element ndmas ftype=integer pytype=int
        
        
        Defined at global.fpp line 351
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ndmas(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ndmas = self._arrays[array_handle]
        else:
            ndmas = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ndmas)
            self._arrays[array_handle] = ndmas
        return ndmas
    
    @ndmas.setter
    def ndmas(self, ndmas):
        self.ndmas[...] = ndmas
    
    @property
    def dmg(self):
        """
        Element dmg ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 352
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dmg(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dmg = self._arrays[array_handle]
        else:
            dmg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dmg)
            self._arrays[array_handle] = dmg
        return dmg
    
    @dmg.setter
    def dmg(self, dmg):
        self.dmg[...] = dmg
    
    @property
    def solution(self):
        """
        Element solution ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 353
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__solution(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            solution = self._arrays[array_handle]
        else:
            solution = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__solution)
            self._arrays[array_handle] = solution
        return solution
    
    @solution.setter
    def solution(self, solution):
        self.solution[...] = solution
    
    @property
    def gmreslastsolution(self):
        """
        Element gmreslastsolution ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 354
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gmreslastsolution(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gmreslastsolution = self._arrays[array_handle]
        else:
            gmreslastsolution = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gmreslastsolution)
            self._arrays[array_handle] = gmreslastsolution
        return gmreslastsolution
    
    @gmreslastsolution.setter
    def gmreslastsolution(self, gmreslastsolution):
        self.gmreslastsolution[...] = gmreslastsolution
    
    @property
    def mbpsi(self):
        """
        Element mbpsi ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 356
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__mbpsi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            mbpsi = self._arrays[array_handle]
        else:
            mbpsi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__mbpsi)
            self._arrays[array_handle] = mbpsi
        return mbpsi
    
    @mbpsi.setter
    def mbpsi(self, mbpsi):
        self.mbpsi[...] = mbpsi
    
    @property
    def liluprecond(self):
        """
        Element liluprecond ftype=logical pytype=bool
        
        
        Defined at global.fpp line 359
        
        """
        return _spec.f90wrap_allglobal__get__liluprecond()
    
    @liluprecond.setter
    def liluprecond(self, liluprecond):
        _spec.f90wrap_allglobal__set__liluprecond(liluprecond)
    
    @property
    def beltramiinverse(self):
        """
        Element beltramiinverse ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 360
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__beltramiinverse(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            beltramiinverse = self._arrays[array_handle]
        else:
            beltramiinverse = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__beltramiinverse)
            self._arrays[array_handle] = beltramiinverse
        return beltramiinverse
    
    @beltramiinverse.setter
    def beltramiinverse(self, beltramiinverse):
        self.beltramiinverse[...] = beltramiinverse
    
    @property
    def diotadxup(self):
        """
        Element diotadxup ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 362
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__diotadxup(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            diotadxup = self._arrays[array_handle]
        else:
            diotadxup = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__diotadxup)
            self._arrays[array_handle] = diotadxup
        return diotadxup
    
    @diotadxup.setter
    def diotadxup(self, diotadxup):
        self.diotadxup[...] = diotadxup
    
    @property
    def ditgpdxtp(self):
        """
        Element ditgpdxtp ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 363
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ditgpdxtp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ditgpdxtp = self._arrays[array_handle]
        else:
            ditgpdxtp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ditgpdxtp)
            self._arrays[array_handle] = ditgpdxtp
        return ditgpdxtp
    
    @ditgpdxtp.setter
    def ditgpdxtp(self, ditgpdxtp):
        self.ditgpdxtp[...] = ditgpdxtp
    
    @property
    def glambda(self):
        """
        Element glambda ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 364
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__glambda(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            glambda = self._arrays[array_handle]
        else:
            glambda = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__glambda)
            self._arrays[array_handle] = glambda
        return glambda
    
    @glambda.setter
    def glambda(self, glambda):
        self.glambda[...] = glambda
    
    @property
    def lmns(self):
        """
        Element lmns ftype=integer               pytype=int
        
        
        Defined at global.fpp line 365
        
        """
        return _spec.f90wrap_allglobal__get__lmns()
    
    @lmns.setter
    def lmns(self, lmns):
        _spec.f90wrap_allglobal__set__lmns(lmns)
    
    @property
    def bemn(self):
        """
        Element bemn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 371
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bemn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bemn = self._arrays[array_handle]
        else:
            bemn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bemn)
            self._arrays[array_handle] = bemn
        return bemn
    
    @bemn.setter
    def bemn(self, bemn):
        self.bemn[...] = bemn
    
    @property
    def iomn(self):
        """
        Element iomn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 371
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iomn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iomn = self._arrays[array_handle]
        else:
            iomn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iomn)
            self._arrays[array_handle] = iomn
        return iomn
    
    @iomn.setter
    def iomn(self, iomn):
        self.iomn[...] = iomn
    
    @property
    def somn(self):
        """
        Element somn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 371
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__somn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            somn = self._arrays[array_handle]
        else:
            somn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__somn)
            self._arrays[array_handle] = somn
        return somn
    
    @somn.setter
    def somn(self, somn):
        self.somn[...] = somn
    
    @property
    def pomn(self):
        """
        Element pomn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 371
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__pomn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pomn = self._arrays[array_handle]
        else:
            pomn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__pomn)
            self._arrays[array_handle] = pomn
        return pomn
    
    @pomn.setter
    def pomn(self, pomn):
        self.pomn[...] = pomn
    
    @property
    def bomn(self):
        """
        Element bomn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 372
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bomn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bomn = self._arrays[array_handle]
        else:
            bomn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bomn)
            self._arrays[array_handle] = bomn
        return bomn
    
    @bomn.setter
    def bomn(self, bomn):
        self.bomn[...] = bomn
    
    @property
    def iemn(self):
        """
        Element iemn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 372
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iemn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iemn = self._arrays[array_handle]
        else:
            iemn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iemn)
            self._arrays[array_handle] = iemn
        return iemn
    
    @iemn.setter
    def iemn(self, iemn):
        self.iemn[...] = iemn
    
    @property
    def semn(self):
        """
        Element semn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 372
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__semn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            semn = self._arrays[array_handle]
        else:
            semn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__semn)
            self._arrays[array_handle] = semn
        return semn
    
    @semn.setter
    def semn(self, semn):
        self.semn[...] = semn
    
    @property
    def pemn(self):
        """
        Element pemn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 372
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__pemn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pemn = self._arrays[array_handle]
        else:
            pemn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__pemn)
            self._arrays[array_handle] = pemn
        return pemn
    
    @pemn.setter
    def pemn(self, pemn):
        self.pemn[...] = pemn
    
    @property
    def bbe(self):
        """
        Element bbe ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 373
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bbe(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bbe = self._arrays[array_handle]
        else:
            bbe = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bbe)
            self._arrays[array_handle] = bbe
        return bbe
    
    @bbe.setter
    def bbe(self, bbe):
        self.bbe[...] = bbe
    
    @property
    def iio(self):
        """
        Element iio ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 373
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iio(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iio = self._arrays[array_handle]
        else:
            iio = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iio)
            self._arrays[array_handle] = iio
        return iio
    
    @iio.setter
    def iio(self, iio):
        self.iio[...] = iio
    
    @property
    def bbo(self):
        """
        Element bbo ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 373
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bbo(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bbo = self._arrays[array_handle]
        else:
            bbo = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bbo)
            self._arrays[array_handle] = bbo
        return bbo
    
    @bbo.setter
    def bbo(self, bbo):
        self.bbo[...] = bbo
    
    @property
    def iie(self):
        """
        Element iie ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 373
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iie(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iie = self._arrays[array_handle]
        else:
            iie = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iie)
            self._arrays[array_handle] = iie
        return iie
    
    @iie.setter
    def iie(self, iie):
        self.iie[...] = iie
    
    @property
    def btemn(self):
        """
        Element btemn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 379
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__btemn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            btemn = self._arrays[array_handle]
        else:
            btemn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__btemn)
            self._arrays[array_handle] = btemn
        return btemn
    
    @btemn.setter
    def btemn(self, btemn):
        self.btemn[...] = btemn
    
    @property
    def bzemn(self):
        """
        Element bzemn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 379
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bzemn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bzemn = self._arrays[array_handle]
        else:
            bzemn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bzemn)
            self._arrays[array_handle] = bzemn
        return bzemn
    
    @bzemn.setter
    def bzemn(self, bzemn):
        self.bzemn[...] = bzemn
    
    @property
    def btomn(self):
        """
        Element btomn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 379
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__btomn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            btomn = self._arrays[array_handle]
        else:
            btomn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__btomn)
            self._arrays[array_handle] = btomn
        return btomn
    
    @btomn.setter
    def btomn(self, btomn):
        self.btomn[...] = btomn
    
    @property
    def bzomn(self):
        """
        Element bzomn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 379
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bzomn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bzomn = self._arrays[array_handle]
        else:
            bzomn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bzomn)
            self._arrays[array_handle] = bzomn
        return bzomn
    
    @bzomn.setter
    def bzomn(self, bzomn):
        self.bzomn[...] = bzomn
    
    @property
    def bloweremn(self):
        """
        Element bloweremn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 385
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bloweremn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bloweremn = self._arrays[array_handle]
        else:
            bloweremn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bloweremn)
            self._arrays[array_handle] = bloweremn
        return bloweremn
    
    @bloweremn.setter
    def bloweremn(self, bloweremn):
        self.bloweremn[...] = bloweremn
    
    @property
    def bloweromn(self):
        """
        Element bloweromn ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 385
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__bloweromn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bloweromn = self._arrays[array_handle]
        else:
            bloweromn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__bloweromn)
            self._arrays[array_handle] = bloweromn
        return bloweromn
    
    @bloweromn.setter
    def bloweromn(self, bloweromn):
        self.bloweromn[...] = bloweromn
    
    @property
    def lgdof(self):
        """
        Element lgdof ftype=integer               pytype=int
        
        
        Defined at global.fpp line 391
        
        """
        return _spec.f90wrap_allglobal__get__lgdof()
    
    @lgdof.setter
    def lgdof(self, lgdof):
        _spec.f90wrap_allglobal__set__lgdof(lgdof)
    
    @property
    def ngdof(self):
        """
        Element ngdof ftype=integer               pytype=int
        
        
        Defined at global.fpp line 392
        
        """
        return _spec.f90wrap_allglobal__get__ngdof()
    
    @ngdof.setter
    def ngdof(self, ngdof):
        _spec.f90wrap_allglobal__set__ngdof(ngdof)
    
    @property
    def dbbdrz(self):
        """
        Element dbbdrz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 400
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dbbdrz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dbbdrz = self._arrays[array_handle]
        else:
            dbbdrz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dbbdrz)
            self._arrays[array_handle] = dbbdrz
        return dbbdrz
    
    @dbbdrz.setter
    def dbbdrz(self, dbbdrz):
        self.dbbdrz[...] = dbbdrz
    
    @property
    def diidrz(self):
        """
        Element diidrz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 401
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__diidrz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            diidrz = self._arrays[array_handle]
        else:
            diidrz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__diidrz)
            self._arrays[array_handle] = diidrz
        return diidrz
    
    @diidrz.setter
    def diidrz(self, diidrz):
        self.diidrz[...] = diidrz
    
    @property
    def dffdrz(self):
        """
        Element dffdrz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 402
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dffdrz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dffdrz = self._arrays[array_handle]
        else:
            dffdrz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dffdrz)
            self._arrays[array_handle] = dffdrz
        return dffdrz
    
    @dffdrz.setter
    def dffdrz(self, dffdrz):
        self.dffdrz[...] = dffdrz
    
    @property
    def dbbdmp(self):
        """
        Element dbbdmp ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 403
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dbbdmp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dbbdmp = self._arrays[array_handle]
        else:
            dbbdmp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dbbdmp)
            self._arrays[array_handle] = dbbdmp
        return dbbdmp
    
    @dbbdmp.setter
    def dbbdmp(self, dbbdmp):
        self.dbbdmp[...] = dbbdmp
    
    @property
    def dmupfdx(self):
        """
        Element dmupfdx ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 445
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dmupfdx(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dmupfdx = self._arrays[array_handle]
        else:
            dmupfdx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dmupfdx)
            self._arrays[array_handle] = dmupfdx
        return dmupfdx
    
    @dmupfdx.setter
    def dmupfdx(self, dmupfdx):
        self.dmupfdx[...] = dmupfdx
    
    @property
    def lhessianallocated(self):
        """
        Element lhessianallocated ftype=logical pytype=bool
        
        
        Defined at global.fpp line 453
        
        """
        return _spec.f90wrap_allglobal__get__lhessianallocated()
    
    @lhessianallocated.setter
    def lhessianallocated(self, lhessianallocated):
        _spec.f90wrap_allglobal__set__lhessianallocated(lhessianallocated)
    
    @property
    def hessian(self):
        """
        Element hessian ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 454
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__hessian(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            hessian = self._arrays[array_handle]
        else:
            hessian = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__hessian)
            self._arrays[array_handle] = hessian
        return hessian
    
    @hessian.setter
    def hessian(self, hessian):
        self.hessian[...] = hessian
    
    @property
    def dessian(self):
        """
        Element dessian ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 455
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dessian(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dessian = self._arrays[array_handle]
        else:
            dessian = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dessian)
            self._arrays[array_handle] = dessian
        return dessian
    
    @dessian.setter
    def dessian(self, dessian):
        self.dessian[...] = dessian
    
    @property
    def cosi(self):
        """
        Element cosi ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 462
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__cosi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cosi = self._arrays[array_handle]
        else:
            cosi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__cosi)
            self._arrays[array_handle] = cosi
        return cosi
    
    @cosi.setter
    def cosi(self, cosi):
        self.cosi[...] = cosi
    
    @property
    def sini(self):
        """
        Element sini ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 462
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__sini(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sini = self._arrays[array_handle]
        else:
            sini = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__sini)
            self._arrays[array_handle] = sini
        return sini
    
    @sini.setter
    def sini(self, sini):
        self.sini[...] = sini
    
    @property
    def gteta(self):
        """
        Element gteta ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 462
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gteta(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gteta = self._arrays[array_handle]
        else:
            gteta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gteta)
            self._arrays[array_handle] = gteta
        return gteta
    
    @gteta.setter
    def gteta(self, gteta):
        self.gteta[...] = gteta
    
    @property
    def gzeta(self):
        """
        Element gzeta ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 462
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gzeta(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gzeta = self._arrays[array_handle]
        else:
            gzeta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gzeta)
            self._arrays[array_handle] = gzeta
        return gzeta
    
    @gzeta.setter
    def gzeta(self, gzeta):
        self.gzeta[...] = gzeta
    
    @property
    def ajk(self):
        """
        Element ajk ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 463
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__ajk(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ajk = self._arrays[array_handle]
        else:
            ajk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__ajk)
            self._arrays[array_handle] = ajk
        return ajk
    
    @ajk.setter
    def ajk(self, ajk):
        self.ajk[...] = ajk
    
    @property
    def dradr(self):
        """
        Element dradr ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 464
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dradr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dradr = self._arrays[array_handle]
        else:
            dradr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dradr)
            self._arrays[array_handle] = dradr
        return dradr
    
    @dradr.setter
    def dradr(self, dradr):
        self.dradr[...] = dradr
    
    @property
    def dradz(self):
        """
        Element dradz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 464
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dradz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dradz = self._arrays[array_handle]
        else:
            dradz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dradz)
            self._arrays[array_handle] = dradz
        return dradz
    
    @dradz.setter
    def dradz(self, dradz):
        self.dradz[...] = dradz
    
    @property
    def dzadr(self):
        """
        Element dzadr ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 464
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzadr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzadr = self._arrays[array_handle]
        else:
            dzadr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzadr)
            self._arrays[array_handle] = dzadr
        return dzadr
    
    @dzadr.setter
    def dzadr(self, dzadr):
        self.dzadr[...] = dzadr
    
    @property
    def dzadz(self):
        """
        Element dzadz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 464
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzadz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzadz = self._arrays[array_handle]
        else:
            dzadz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzadz)
            self._arrays[array_handle] = dzadz
        return dzadz
    
    @dzadz.setter
    def dzadz(self, dzadz):
        self.dzadz[...] = dzadz
    
    @property
    def drodr(self):
        """
        Element drodr ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 465
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__drodr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            drodr = self._arrays[array_handle]
        else:
            drodr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__drodr)
            self._arrays[array_handle] = drodr
        return drodr
    
    @drodr.setter
    def drodr(self, drodr):
        self.drodr[...] = drodr
    
    @property
    def drodz(self):
        """
        Element drodz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 465
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__drodz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            drodz = self._arrays[array_handle]
        else:
            drodz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__drodz)
            self._arrays[array_handle] = drodz
        return drodz
    
    @drodz.setter
    def drodz(self, drodz):
        self.drodz[...] = drodz
    
    @property
    def dzodr(self):
        """
        Element dzodr ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 465
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzodr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzodr = self._arrays[array_handle]
        else:
            dzodr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzodr)
            self._arrays[array_handle] = dzodr
        return dzodr
    
    @dzodr.setter
    def dzodr(self, dzodr):
        self.dzodr[...] = dzodr
    
    @property
    def dzodz(self):
        """
        Element dzodz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 465
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dzodz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dzodz = self._arrays[array_handle]
        else:
            dzodz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dzodz)
            self._arrays[array_handle] = dzodz
        return dzodz
    
    @dzodz.setter
    def dzodz(self, dzodz):
        self.dzodz[...] = dzodz
    
    @property
    def djkp(self):
        """
        Element djkp ftype=integer pytype=int
        
        
        Defined at global.fpp line 466
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__djkp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            djkp = self._arrays[array_handle]
        else:
            djkp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__djkp)
            self._arrays[array_handle] = djkp
        return djkp
    
    @djkp.setter
    def djkp(self, djkp):
        self.djkp[...] = djkp
    
    @property
    def djkm(self):
        """
        Element djkm ftype=integer pytype=int
        
        
        Defined at global.fpp line 466
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__djkm(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            djkm = self._arrays[array_handle]
        else:
            djkm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__djkm)
            self._arrays[array_handle] = djkm
        return djkm
    
    @djkm.setter
    def djkm(self, djkm):
        self.djkm[...] = djkm
    
    @property
    def lbbintegral(self):
        """
        Element lbbintegral ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 498
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__lbbintegral(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lbbintegral = self._arrays[array_handle]
        else:
            lbbintegral = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__lbbintegral)
            self._arrays[array_handle] = lbbintegral
        return lbbintegral
    
    @lbbintegral.setter
    def lbbintegral(self, lbbintegral):
        self.lbbintegral[...] = lbbintegral
    
    @property
    def labintegral(self):
        """
        Element labintegral ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 499
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__labintegral(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            labintegral = self._arrays[array_handle]
        else:
            labintegral = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__labintegral)
            self._arrays[array_handle] = labintegral
        return labintegral
    
    @labintegral.setter
    def labintegral(self, labintegral):
        self.labintegral[...] = labintegral
    
    @property
    def vvolume(self):
        """
        Element vvolume ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 503
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__vvolume(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            vvolume = self._arrays[array_handle]
        else:
            vvolume = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__vvolume)
            self._arrays[array_handle] = vvolume
        return vvolume
    
    @vvolume.setter
    def vvolume(self, vvolume):
        self.vvolume[...] = vvolume
    
    @property
    def dvolume(self):
        """
        Element dvolume ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 504
        
        """
        return _spec.f90wrap_allglobal__get__dvolume()
    
    @dvolume.setter
    def dvolume(self, dvolume):
        _spec.f90wrap_allglobal__set__dvolume(dvolume)
    
    @property
    def ivol(self):
        """
        Element ivol ftype=integer               pytype=int
        
        
        Defined at global.fpp line 507
        
        """
        return _spec.f90wrap_allglobal__get__ivol()
    
    @ivol.setter
    def ivol(self, ivol):
        _spec.f90wrap_allglobal__set__ivol(ivol)
    
    @property
    def gbzeta(self):
        """
        Element gbzeta ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 508
        
        """
        return _spec.f90wrap_allglobal__get__gbzeta()
    
    @gbzeta.setter
    def gbzeta(self, gbzeta):
        _spec.f90wrap_allglobal__set__gbzeta(gbzeta)
    
    @property
    def iquad(self):
        """
        Element iquad ftype=integer pytype=int
        
        
        Defined at global.fpp line 509
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__iquad(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iquad = self._arrays[array_handle]
        else:
            iquad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__iquad)
            self._arrays[array_handle] = iquad
        return iquad
    
    @iquad.setter
    def iquad(self, iquad):
        self.iquad[...] = iquad
    
    @property
    def gaussianweight(self):
        """
        Element gaussianweight ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 510
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gaussianweight(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gaussianweight = self._arrays[array_handle]
        else:
            gaussianweight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gaussianweight)
            self._arrays[array_handle] = gaussianweight
        return gaussianweight
    
    @gaussianweight.setter
    def gaussianweight(self, gaussianweight):
        self.gaussianweight[...] = gaussianweight
    
    @property
    def gaussianabscissae(self):
        """
        Element gaussianabscissae ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 510
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__gaussianabscissae(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gaussianabscissae = self._arrays[array_handle]
        else:
            gaussianabscissae = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__gaussianabscissae)
            self._arrays[array_handle] = gaussianabscissae
        return gaussianabscissae
    
    @gaussianabscissae.setter
    def gaussianabscissae(self, gaussianabscissae):
        self.gaussianabscissae[...] = gaussianabscissae
    
    @property
    def lblinear(self):
        """
        Element lblinear ftype=logical pytype=bool
        
        
        Defined at global.fpp line 511
        
        """
        return _spec.f90wrap_allglobal__get__lblinear()
    
    @lblinear.setter
    def lblinear(self, lblinear):
        _spec.f90wrap_allglobal__set__lblinear(lblinear)
    
    @property
    def lbnewton(self):
        """
        Element lbnewton ftype=logical pytype=bool
        
        
        Defined at global.fpp line 511
        
        """
        return _spec.f90wrap_allglobal__get__lbnewton()
    
    @lbnewton.setter
    def lbnewton(self, lbnewton):
        _spec.f90wrap_allglobal__set__lbnewton(lbnewton)
    
    @property
    def lbsequad(self):
        """
        Element lbsequad ftype=logical pytype=bool
        
        
        Defined at global.fpp line 511
        
        """
        return _spec.f90wrap_allglobal__get__lbsequad()
    
    @lbsequad.setter
    def lbsequad(self, lbsequad):
        _spec.f90wrap_allglobal__set__lbsequad(lbsequad)
    
    @property
    def orzp(self):
        """
        Element orzp ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 512
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__orzp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            orzp = self._arrays[array_handle]
        else:
            orzp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__orzp)
            self._arrays[array_handle] = orzp
        return orzp
    
    @orzp.setter
    def orzp(self, orzp):
        self.orzp[...] = orzp
    
    @property
    def globaljk(self):
        """
        Element globaljk ftype=integer               pytype=int
        
        
        Defined at global.fpp line 520
        
        """
        return _spec.f90wrap_allglobal__get__globaljk()
    
    @globaljk.setter
    def globaljk(self, globaljk):
        _spec.f90wrap_allglobal__set__globaljk(globaljk)
    
    @property
    def dxyz(self):
        """
        Element dxyz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 521
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__dxyz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dxyz = self._arrays[array_handle]
        else:
            dxyz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__dxyz)
            self._arrays[array_handle] = dxyz
        return dxyz
    
    @dxyz.setter
    def dxyz(self, dxyz):
        self.dxyz[...] = dxyz
    
    @property
    def nxyz(self):
        """
        Element nxyz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 522
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__nxyz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            nxyz = self._arrays[array_handle]
        else:
            nxyz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__nxyz)
            self._arrays[array_handle] = nxyz
        return nxyz
    
    @nxyz.setter
    def nxyz(self, nxyz):
        self.nxyz[...] = nxyz
    
    @property
    def jxyz(self):
        """
        Element jxyz ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 523
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__jxyz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jxyz = self._arrays[array_handle]
        else:
            jxyz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__jxyz)
            self._arrays[array_handle] = jxyz
        return jxyz
    
    @jxyz.setter
    def jxyz(self, jxyz):
        self.jxyz[...] = jxyz
    
    @property
    def tetazeta(self):
        """
        Element tetazeta ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 524
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_allglobal__array__tetazeta(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tetazeta = self._arrays[array_handle]
        else:
            tetazeta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_allglobal__array__tetazeta)
            self._arrays[array_handle] = tetazeta
        return tetazeta
    
    @tetazeta.setter
    def tetazeta(self, tetazeta):
        self.tetazeta[...] = tetazeta
    
    @property
    def virtualcasingfactor(self):
        """
        Element virtualcasingfactor ftype=real(8) pytype=float
        
        
        Defined at global.fpp line 526
        
        """
        return _spec.f90wrap_allglobal__get__virtualcasingfactor()
    
    @virtualcasingfactor.setter
    def virtualcasingfactor(self, virtualcasingfactor):
        _spec.f90wrap_allglobal__set__virtualcasingfactor(virtualcasingfactor)
    
    @property
    def iberror(self):
        """
        Element iberror ftype=integer               pytype=int
        
        
        Defined at global.fpp line 527
        
        """
        return _spec.f90wrap_allglobal__get__iberror()
    
    @iberror.setter
    def iberror(self, iberror):
        _spec.f90wrap_allglobal__set__iberror(iberror)
    
    @property
    def nfreeboundaryiterations(self):
        """
        Element nfreeboundaryiterations ftype=integer               pytype=int
        
        
        Defined at global.fpp line 528
        
        """
        return _spec.f90wrap_allglobal__get__nfreeboundaryiterations()
    
    @nfreeboundaryiterations.setter
    def nfreeboundaryiterations(self, nfreeboundaryiterations):
        _spec.f90wrap_allglobal__set__nfreeboundaryiterations(nfreeboundaryiterations)
    
    @property
    def node(self):
        """
        Element node ftype=integer pytype=int
        
        
        Defined at global.fpp line 530
        
        """
        return _spec.f90wrap_allglobal__get__node()
    
    @property
    def first_free_bound(self):
        """
        Element first_free_bound ftype=logical pytype=bool
        
        
        Defined at global.fpp line 532
        
        """
        return _spec.f90wrap_allglobal__get__first_free_bound()
    
    @first_free_bound.setter
    def first_free_bound(self, first_free_bound):
        _spec.f90wrap_allglobal__set__first_free_bound(first_free_bound)
    
    def __str__(self):
        ret = ['<allglobal>{\n']
        ret.append('    myid : ')
        ret.append(repr(self.myid))
        ret.append(',\n    ncpu : ')
        ret.append(repr(self.ncpu))
        ret.append(',\n    mpi_comm_spec : ')
        ret.append(repr(self.mpi_comm_spec))
        ret.append(',\n    ismyvolumevalue : ')
        ret.append(repr(self.ismyvolumevalue))
        ret.append(',\n    cpus : ')
        ret.append(repr(self.cpus))
        ret.append(',\n    skip_write : ')
        ret.append(repr(self.skip_write))
        ret.append(',\n    pi2nfp : ')
        ret.append(repr(self.pi2nfp))
        ret.append(',\n    pi2pi2nfp : ')
        ret.append(repr(self.pi2pi2nfp))
        ret.append(',\n    pi2pi2nfphalf : ')
        ret.append(repr(self.pi2pi2nfphalf))
        ret.append(',\n    pi2pi2nfpquart : ')
        ret.append(repr(self.pi2pi2nfpquart))
        ret.append(',\n    ext : ')
        ret.append(repr(self.ext))
        ret.append(',\n    forceerr : ')
        ret.append(repr(self.forceerr))
        ret.append(',\n    energy : ')
        ret.append(repr(self.energy))
        ret.append(',\n    ipdt : ')
        ret.append(repr(self.ipdt))
        ret.append(',\n    ipdtdpf : ')
        ret.append(repr(self.ipdtdpf))
        ret.append(',\n    mvol : ')
        ret.append(repr(self.mvol))
        ret.append(',\n    yesstellsym : ')
        ret.append(repr(self.yesstellsym))
        ret.append(',\n    notstellsym : ')
        ret.append(repr(self.notstellsym))
        ret.append(',\n    yesmatrixfree : ')
        ret.append(repr(self.yesmatrixfree))
        ret.append(',\n    notmatrixfree : ')
        ret.append(repr(self.notmatrixfree))
        ret.append(',\n    cheby : ')
        ret.append(repr(self.cheby))
        ret.append(',\n    zernike : ')
        ret.append(repr(self.zernike))
        ret.append(',\n    tt : ')
        ret.append(repr(self.tt))
        ret.append(',\n    rtt : ')
        ret.append(repr(self.rtt))
        ret.append(',\n    rtm : ')
        ret.append(repr(self.rtm))
        ret.append(',\n    zernikedof : ')
        ret.append(repr(self.zernikedof))
        ret.append(',\n    imagneticok : ')
        ret.append(repr(self.imagneticok))
        ret.append(',\n    iconstraintok : ')
        ret.append(repr(self.iconstraintok))
        ret.append(',\n    beltramierror : ')
        ret.append(repr(self.beltramierror))
        ret.append(',\n    mn : ')
        ret.append(repr(self.mn))
        ret.append(',\n    im : ')
        ret.append(repr(self.im))
        ret.append(',\n    in_ : ')
        ret.append(repr(self.in_))
        ret.append(',\n    halfmm : ')
        ret.append(repr(self.halfmm))
        ret.append(',\n    regumm : ')
        ret.append(repr(self.regumm))
        ret.append(',\n    rscale : ')
        ret.append(repr(self.rscale))
        ret.append(',\n    psifactor : ')
        ret.append(repr(self.psifactor))
        ret.append(',\n    inifactor : ')
        ret.append(repr(self.inifactor))
        ret.append(',\n    bbweight : ')
        ret.append(repr(self.bbweight))
        ret.append(',\n    mmpp : ')
        ret.append(repr(self.mmpp))
        ret.append(',\n    mne : ')
        ret.append(repr(self.mne))
        ret.append(',\n    ime : ')
        ret.append(repr(self.ime))
        ret.append(',\n    ine : ')
        ret.append(repr(self.ine))
        ret.append(',\n    mns : ')
        ret.append(repr(self.mns))
        ret.append(',\n    ims : ')
        ret.append(repr(self.ims))
        ret.append(',\n    ins : ')
        ret.append(repr(self.ins))
        ret.append(',\n    lmpol : ')
        ret.append(repr(self.lmpol))
        ret.append(',\n    lntor : ')
        ret.append(repr(self.lntor))
        ret.append(',\n    smpol : ')
        ret.append(repr(self.smpol))
        ret.append(',\n    sntor : ')
        ret.append(repr(self.sntor))
        ret.append(',\n    xoffset : ')
        ret.append(repr(self.xoffset))
        ret.append(',\n    irbc : ')
        ret.append(repr(self.irbc))
        ret.append(',\n    izbs : ')
        ret.append(repr(self.izbs))
        ret.append(',\n    irbs : ')
        ret.append(repr(self.irbs))
        ret.append(',\n    izbc : ')
        ret.append(repr(self.izbc))
        ret.append(',\n    drbc : ')
        ret.append(repr(self.drbc))
        ret.append(',\n    dzbs : ')
        ret.append(repr(self.dzbs))
        ret.append(',\n    drbs : ')
        ret.append(repr(self.drbs))
        ret.append(',\n    dzbc : ')
        ret.append(repr(self.dzbc))
        ret.append(',\n    irij : ')
        ret.append(repr(self.irij))
        ret.append(',\n    izij : ')
        ret.append(repr(self.izij))
        ret.append(',\n    drij : ')
        ret.append(repr(self.drij))
        ret.append(',\n    dzij : ')
        ret.append(repr(self.dzij))
        ret.append(',\n    trij : ')
        ret.append(repr(self.trij))
        ret.append(',\n    tzij : ')
        ret.append(repr(self.tzij))
        ret.append(',\n    ivns : ')
        ret.append(repr(self.ivns))
        ret.append(',\n    ibns : ')
        ret.append(repr(self.ibns))
        ret.append(',\n    ivnc : ')
        ret.append(repr(self.ivnc))
        ret.append(',\n    ibnc : ')
        ret.append(repr(self.ibnc))
        ret.append(',\n    lrbc : ')
        ret.append(repr(self.lrbc))
        ret.append(',\n    lzbs : ')
        ret.append(repr(self.lzbs))
        ret.append(',\n    lrbs : ')
        ret.append(repr(self.lrbs))
        ret.append(',\n    lzbc : ')
        ret.append(repr(self.lzbc))
        ret.append(',\n    num_modes : ')
        ret.append(repr(self.num_modes))
        ret.append(',\n    mmrzrz : ')
        ret.append(repr(self.mmrzrz))
        ret.append(',\n    nnrzrz : ')
        ret.append(repr(self.nnrzrz))
        ret.append(',\n    allrzrz : ')
        ret.append(repr(self.allrzrz))
        ret.append(',\n    nt : ')
        ret.append(repr(self.nt))
        ret.append(',\n    nz : ')
        ret.append(repr(self.nz))
        ret.append(',\n    ntz : ')
        ret.append(repr(self.ntz))
        ret.append(',\n    hnt : ')
        ret.append(repr(self.hnt))
        ret.append(',\n    hnz : ')
        ret.append(repr(self.hnz))
        ret.append(',\n    sontz : ')
        ret.append(repr(self.sontz))
        ret.append(',\n    rij : ')
        ret.append(repr(self.rij))
        ret.append(',\n    zij : ')
        ret.append(repr(self.zij))
        ret.append(',\n    xij : ')
        ret.append(repr(self.xij))
        ret.append(',\n    yij : ')
        ret.append(repr(self.yij))
        ret.append(',\n    sg : ')
        ret.append(repr(self.sg))
        ret.append(',\n    guvij : ')
        ret.append(repr(self.guvij))
        ret.append(',\n    gvuij : ')
        ret.append(repr(self.gvuij))
        ret.append(',\n    guvijsave : ')
        ret.append(repr(self.guvijsave))
        ret.append(',\n    ki : ')
        ret.append(repr(self.ki))
        ret.append(',\n    kijs : ')
        ret.append(repr(self.kijs))
        ret.append(',\n    kija : ')
        ret.append(repr(self.kija))
        ret.append(',\n    iotakkii : ')
        ret.append(repr(self.iotakkii))
        ret.append(',\n    iotaksub : ')
        ret.append(repr(self.iotaksub))
        ret.append(',\n    iotakadd : ')
        ret.append(repr(self.iotakadd))
        ret.append(',\n    iotaksgn : ')
        ret.append(repr(self.iotaksgn))
        ret.append(',\n    efmn : ')
        ret.append(repr(self.efmn))
        ret.append(',\n    ofmn : ')
        ret.append(repr(self.ofmn))
        ret.append(',\n    cfmn : ')
        ret.append(repr(self.cfmn))
        ret.append(',\n    sfmn : ')
        ret.append(repr(self.sfmn))
        ret.append(',\n    evmn : ')
        ret.append(repr(self.evmn))
        ret.append(',\n    odmn : ')
        ret.append(repr(self.odmn))
        ret.append(',\n    comn : ')
        ret.append(repr(self.comn))
        ret.append(',\n    simn : ')
        ret.append(repr(self.simn))
        ret.append(',\n    ijreal : ')
        ret.append(repr(self.ijreal))
        ret.append(',\n    ijimag : ')
        ret.append(repr(self.ijimag))
        ret.append(',\n    jireal : ')
        ret.append(repr(self.jireal))
        ret.append(',\n    jiimag : ')
        ret.append(repr(self.jiimag))
        ret.append(',\n    jkreal : ')
        ret.append(repr(self.jkreal))
        ret.append(',\n    jkimag : ')
        ret.append(repr(self.jkimag))
        ret.append(',\n    kjreal : ')
        ret.append(repr(self.kjreal))
        ret.append(',\n    kjimag : ')
        ret.append(repr(self.kjimag))
        ret.append(',\n    bsupumn : ')
        ret.append(repr(self.bsupumn))
        ret.append(',\n    bsupvmn : ')
        ret.append(repr(self.bsupvmn))
        ret.append(',\n    goomne : ')
        ret.append(repr(self.goomne))
        ret.append(',\n    goomno : ')
        ret.append(repr(self.goomno))
        ret.append(',\n    gssmne : ')
        ret.append(repr(self.gssmne))
        ret.append(',\n    gssmno : ')
        ret.append(repr(self.gssmno))
        ret.append(',\n    gstmne : ')
        ret.append(repr(self.gstmne))
        ret.append(',\n    gstmno : ')
        ret.append(repr(self.gstmno))
        ret.append(',\n    gszmne : ')
        ret.append(repr(self.gszmne))
        ret.append(',\n    gszmno : ')
        ret.append(repr(self.gszmno))
        ret.append(',\n    gttmne : ')
        ret.append(repr(self.gttmne))
        ret.append(',\n    gttmno : ')
        ret.append(repr(self.gttmno))
        ret.append(',\n    gtzmne : ')
        ret.append(repr(self.gtzmne))
        ret.append(',\n    gtzmno : ')
        ret.append(repr(self.gtzmno))
        ret.append(',\n    gzzmne : ')
        ret.append(repr(self.gzzmne))
        ret.append(',\n    gzzmno : ')
        ret.append(repr(self.gzzmno))
        ret.append(',\n    dtoocc : ')
        ret.append(repr(self.dtoocc))
        ret.append(',\n    dtoocs : ')
        ret.append(repr(self.dtoocs))
        ret.append(',\n    dtoosc : ')
        ret.append(repr(self.dtoosc))
        ret.append(',\n    dtooss : ')
        ret.append(repr(self.dtooss))
        ret.append(',\n    ttsscc : ')
        ret.append(repr(self.ttsscc))
        ret.append(',\n    ttsscs : ')
        ret.append(repr(self.ttsscs))
        ret.append(',\n    ttsssc : ')
        ret.append(repr(self.ttsssc))
        ret.append(',\n    ttssss : ')
        ret.append(repr(self.ttssss))
        ret.append(',\n    tdstcc : ')
        ret.append(repr(self.tdstcc))
        ret.append(',\n    tdstcs : ')
        ret.append(repr(self.tdstcs))
        ret.append(',\n    tdstsc : ')
        ret.append(repr(self.tdstsc))
        ret.append(',\n    tdstss : ')
        ret.append(repr(self.tdstss))
        ret.append(',\n    tdszcc : ')
        ret.append(repr(self.tdszcc))
        ret.append(',\n    tdszcs : ')
        ret.append(repr(self.tdszcs))
        ret.append(',\n    tdszsc : ')
        ret.append(repr(self.tdszsc))
        ret.append(',\n    tdszss : ')
        ret.append(repr(self.tdszss))
        ret.append(',\n    ddttcc : ')
        ret.append(repr(self.ddttcc))
        ret.append(',\n    ddttcs : ')
        ret.append(repr(self.ddttcs))
        ret.append(',\n    ddttsc : ')
        ret.append(repr(self.ddttsc))
        ret.append(',\n    ddttss : ')
        ret.append(repr(self.ddttss))
        ret.append(',\n    ddtzcc : ')
        ret.append(repr(self.ddtzcc))
        ret.append(',\n    ddtzcs : ')
        ret.append(repr(self.ddtzcs))
        ret.append(',\n    ddtzsc : ')
        ret.append(repr(self.ddtzsc))
        ret.append(',\n    ddtzss : ')
        ret.append(repr(self.ddtzss))
        ret.append(',\n    ddzzcc : ')
        ret.append(repr(self.ddzzcc))
        ret.append(',\n    ddzzcs : ')
        ret.append(repr(self.ddzzcs))
        ret.append(',\n    ddzzsc : ')
        ret.append(repr(self.ddzzsc))
        ret.append(',\n    ddzzss : ')
        ret.append(repr(self.ddzzss))
        ret.append(',\n    tsc : ')
        ret.append(repr(self.tsc))
        ret.append(',\n    tss : ')
        ret.append(repr(self.tss))
        ret.append(',\n    dtc : ')
        ret.append(repr(self.dtc))
        ret.append(',\n    dts : ')
        ret.append(repr(self.dts))
        ret.append(',\n    dzc : ')
        ret.append(repr(self.dzc))
        ret.append(',\n    dzs : ')
        ret.append(repr(self.dzs))
        ret.append(',\n    ttc : ')
        ret.append(repr(self.ttc))
        ret.append(',\n    tzc : ')
        ret.append(repr(self.tzc))
        ret.append(',\n    tts : ')
        ret.append(repr(self.tts))
        ret.append(',\n    tzs : ')
        ret.append(repr(self.tzs))
        ret.append(',\n    dtflux : ')
        ret.append(repr(self.dtflux))
        ret.append(',\n    dpflux : ')
        ret.append(repr(self.dpflux))
        ret.append(',\n    sweight : ')
        ret.append(repr(self.sweight))
        ret.append(',\n    nadof : ')
        ret.append(repr(self.nadof))
        ret.append(',\n    nfielddof : ')
        ret.append(repr(self.nfielddof))
        ret.append(',\n    lma : ')
        ret.append(repr(self.lma))
        ret.append(',\n    lmb : ')
        ret.append(repr(self.lmb))
        ret.append(',\n    lmc : ')
        ret.append(repr(self.lmc))
        ret.append(',\n    lmd : ')
        ret.append(repr(self.lmd))
        ret.append(',\n    lme : ')
        ret.append(repr(self.lme))
        ret.append(',\n    lmf : ')
        ret.append(repr(self.lmf))
        ret.append(',\n    lmg : ')
        ret.append(repr(self.lmg))
        ret.append(',\n    lmh : ')
        ret.append(repr(self.lmh))
        ret.append(',\n    lmavalue : ')
        ret.append(repr(self.lmavalue))
        ret.append(',\n    lmbvalue : ')
        ret.append(repr(self.lmbvalue))
        ret.append(',\n    lmcvalue : ')
        ret.append(repr(self.lmcvalue))
        ret.append(',\n    lmdvalue : ')
        ret.append(repr(self.lmdvalue))
        ret.append(',\n    lmevalue : ')
        ret.append(repr(self.lmevalue))
        ret.append(',\n    lmfvalue : ')
        ret.append(repr(self.lmfvalue))
        ret.append(',\n    lmgvalue : ')
        ret.append(repr(self.lmgvalue))
        ret.append(',\n    lmhvalue : ')
        ret.append(repr(self.lmhvalue))
        ret.append(',\n    fso : ')
        ret.append(repr(self.fso))
        ret.append(',\n    fse : ')
        ret.append(repr(self.fse))
        ret.append(',\n    lcoordinatesingularity : ')
        ret.append(repr(self.lcoordinatesingularity))
        ret.append(',\n    lplasmaregion : ')
        ret.append(repr(self.lplasmaregion))
        ret.append(',\n    lvacuumregion : ')
        ret.append(repr(self.lvacuumregion))
        ret.append(',\n    lsavedguvij : ')
        ret.append(repr(self.lsavedguvij))
        ret.append(',\n    localconstraint : ')
        ret.append(repr(self.localconstraint))
        ret.append(',\n    dma : ')
        ret.append(repr(self.dma))
        ret.append(',\n    dmb : ')
        ret.append(repr(self.dmb))
        ret.append(',\n    dmd : ')
        ret.append(repr(self.dmd))
        ret.append(',\n    dmas : ')
        ret.append(repr(self.dmas))
        ret.append(',\n    dmds : ')
        ret.append(repr(self.dmds))
        ret.append(',\n    idmas : ')
        ret.append(repr(self.idmas))
        ret.append(',\n    jdmas : ')
        ret.append(repr(self.jdmas))
        ret.append(',\n    ndmasmax : ')
        ret.append(repr(self.ndmasmax))
        ret.append(',\n    ndmas : ')
        ret.append(repr(self.ndmas))
        ret.append(',\n    dmg : ')
        ret.append(repr(self.dmg))
        ret.append(',\n    solution : ')
        ret.append(repr(self.solution))
        ret.append(',\n    gmreslastsolution : ')
        ret.append(repr(self.gmreslastsolution))
        ret.append(',\n    mbpsi : ')
        ret.append(repr(self.mbpsi))
        ret.append(',\n    liluprecond : ')
        ret.append(repr(self.liluprecond))
        ret.append(',\n    beltramiinverse : ')
        ret.append(repr(self.beltramiinverse))
        ret.append(',\n    diotadxup : ')
        ret.append(repr(self.diotadxup))
        ret.append(',\n    ditgpdxtp : ')
        ret.append(repr(self.ditgpdxtp))
        ret.append(',\n    glambda : ')
        ret.append(repr(self.glambda))
        ret.append(',\n    lmns : ')
        ret.append(repr(self.lmns))
        ret.append(',\n    bemn : ')
        ret.append(repr(self.bemn))
        ret.append(',\n    iomn : ')
        ret.append(repr(self.iomn))
        ret.append(',\n    somn : ')
        ret.append(repr(self.somn))
        ret.append(',\n    pomn : ')
        ret.append(repr(self.pomn))
        ret.append(',\n    bomn : ')
        ret.append(repr(self.bomn))
        ret.append(',\n    iemn : ')
        ret.append(repr(self.iemn))
        ret.append(',\n    semn : ')
        ret.append(repr(self.semn))
        ret.append(',\n    pemn : ')
        ret.append(repr(self.pemn))
        ret.append(',\n    bbe : ')
        ret.append(repr(self.bbe))
        ret.append(',\n    iio : ')
        ret.append(repr(self.iio))
        ret.append(',\n    bbo : ')
        ret.append(repr(self.bbo))
        ret.append(',\n    iie : ')
        ret.append(repr(self.iie))
        ret.append(',\n    btemn : ')
        ret.append(repr(self.btemn))
        ret.append(',\n    bzemn : ')
        ret.append(repr(self.bzemn))
        ret.append(',\n    btomn : ')
        ret.append(repr(self.btomn))
        ret.append(',\n    bzomn : ')
        ret.append(repr(self.bzomn))
        ret.append(',\n    bloweremn : ')
        ret.append(repr(self.bloweremn))
        ret.append(',\n    bloweromn : ')
        ret.append(repr(self.bloweromn))
        ret.append(',\n    lgdof : ')
        ret.append(repr(self.lgdof))
        ret.append(',\n    ngdof : ')
        ret.append(repr(self.ngdof))
        ret.append(',\n    dbbdrz : ')
        ret.append(repr(self.dbbdrz))
        ret.append(',\n    diidrz : ')
        ret.append(repr(self.diidrz))
        ret.append(',\n    dffdrz : ')
        ret.append(repr(self.dffdrz))
        ret.append(',\n    dbbdmp : ')
        ret.append(repr(self.dbbdmp))
        ret.append(',\n    dmupfdx : ')
        ret.append(repr(self.dmupfdx))
        ret.append(',\n    lhessianallocated : ')
        ret.append(repr(self.lhessianallocated))
        ret.append(',\n    hessian : ')
        ret.append(repr(self.hessian))
        ret.append(',\n    dessian : ')
        ret.append(repr(self.dessian))
        ret.append(',\n    cosi : ')
        ret.append(repr(self.cosi))
        ret.append(',\n    sini : ')
        ret.append(repr(self.sini))
        ret.append(',\n    gteta : ')
        ret.append(repr(self.gteta))
        ret.append(',\n    gzeta : ')
        ret.append(repr(self.gzeta))
        ret.append(',\n    ajk : ')
        ret.append(repr(self.ajk))
        ret.append(',\n    dradr : ')
        ret.append(repr(self.dradr))
        ret.append(',\n    dradz : ')
        ret.append(repr(self.dradz))
        ret.append(',\n    dzadr : ')
        ret.append(repr(self.dzadr))
        ret.append(',\n    dzadz : ')
        ret.append(repr(self.dzadz))
        ret.append(',\n    drodr : ')
        ret.append(repr(self.drodr))
        ret.append(',\n    drodz : ')
        ret.append(repr(self.drodz))
        ret.append(',\n    dzodr : ')
        ret.append(repr(self.dzodr))
        ret.append(',\n    dzodz : ')
        ret.append(repr(self.dzodz))
        ret.append(',\n    djkp : ')
        ret.append(repr(self.djkp))
        ret.append(',\n    djkm : ')
        ret.append(repr(self.djkm))
        ret.append(',\n    lbbintegral : ')
        ret.append(repr(self.lbbintegral))
        ret.append(',\n    labintegral : ')
        ret.append(repr(self.labintegral))
        ret.append(',\n    vvolume : ')
        ret.append(repr(self.vvolume))
        ret.append(',\n    dvolume : ')
        ret.append(repr(self.dvolume))
        ret.append(',\n    ivol : ')
        ret.append(repr(self.ivol))
        ret.append(',\n    gbzeta : ')
        ret.append(repr(self.gbzeta))
        ret.append(',\n    iquad : ')
        ret.append(repr(self.iquad))
        ret.append(',\n    gaussianweight : ')
        ret.append(repr(self.gaussianweight))
        ret.append(',\n    gaussianabscissae : ')
        ret.append(repr(self.gaussianabscissae))
        ret.append(',\n    lblinear : ')
        ret.append(repr(self.lblinear))
        ret.append(',\n    lbnewton : ')
        ret.append(repr(self.lbnewton))
        ret.append(',\n    lbsequad : ')
        ret.append(repr(self.lbsequad))
        ret.append(',\n    orzp : ')
        ret.append(repr(self.orzp))
        ret.append(',\n    globaljk : ')
        ret.append(repr(self.globaljk))
        ret.append(',\n    dxyz : ')
        ret.append(repr(self.dxyz))
        ret.append(',\n    nxyz : ')
        ret.append(repr(self.nxyz))
        ret.append(',\n    jxyz : ')
        ret.append(repr(self.jxyz))
        ret.append(',\n    tetazeta : ')
        ret.append(repr(self.tetazeta))
        ret.append(',\n    virtualcasingfactor : ')
        ret.append(repr(self.virtualcasingfactor))
        ret.append(',\n    iberror : ')
        ret.append(repr(self.iberror))
        ret.append(',\n    nfreeboundaryiterations : ')
        ret.append(repr(self.nfreeboundaryiterations))
        ret.append(',\n    node : ')
        ret.append(repr(self.node))
        ret.append(',\n    first_free_bound : ')
        ret.append(repr(self.first_free_bound))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

allglobal = Allglobal()

class Fftw_Interface(f90wrap.runtime.FortranModule):
    """
    Module fftw_interface
    
    
    Defined at global.fpp lines 2274-2279
    
    """
    @property
    def cplxin(self):
        """
        Element cplxin ftype=complex(c_double_complex) pytype=complex
        
        
        Defined at global.fpp line 2279
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_fftw_interface__array__cplxin(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cplxin = self._arrays[array_handle]
        else:
            cplxin = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_fftw_interface__array__cplxin)
            self._arrays[array_handle] = cplxin
        return cplxin
    
    @cplxin.setter
    def cplxin(self, cplxin):
        self.cplxin[...] = cplxin
    
    @property
    def cplxout(self):
        """
        Element cplxout ftype=complex(c_double_complex) pytype=complex
        
        
        Defined at global.fpp line 2279
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _spec.f90wrap_fftw_interface__array__cplxout(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cplxout = self._arrays[array_handle]
        else:
            cplxout = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _spec.f90wrap_fftw_interface__array__cplxout)
            self._arrays[array_handle] = cplxout
        return cplxout
    
    @cplxout.setter
    def cplxout(self, cplxout):
        self.cplxout[...] = cplxout
    
    def __str__(self):
        ret = ['<fftw_interface>{\n']
        ret.append('    cplxin : ')
        ret.append(repr(self.cplxin))
        ret.append(',\n    cplxout : ')
        ret.append(repr(self.cplxout))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

fftw_interface = Fftw_Interface()

class Intghs_Module(f90wrap.runtime.FortranModule):
    """
    Module intghs_module
    
    
    Defined at intghs.fpp lines 40-49
    
    """
    @f90wrap.runtime.register_class("spec.intghs_workspace")
    class intghs_workspace(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=intghs_workspace)
        
        
        Defined at intghs.fpp lines 41-47
        
        """
        def __init__(self, handle=None):
            """
            self = Intghs_Workspace()
            
            
            Defined at intghs.fpp lines 41-47
            
            
            Returns
            -------
            this : Intghs_Workspace
            	Object to be constructed
            
            
            Automatically generated constructor for intghs_workspace
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _spec.f90wrap_intghs_workspace_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Intghs_Workspace
            
            
            Defined at intghs.fpp lines 41-47
            
            Parameters
            ----------
            this : Intghs_Workspace
            	Object to be destructed
            
            
            Automatically generated destructor for intghs_workspace
            """
            if self._alloc:
                _spec.f90wrap_intghs_workspace_finalise(this=self._handle)
        
        @property
        def efmn(self):
            """
            Element efmn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__efmn(self._handle)
            if array_handle in self._arrays:
                efmn = self._arrays[array_handle]
            else:
                efmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__efmn)
                self._arrays[array_handle] = efmn
            return efmn
        
        @efmn.setter
        def efmn(self, efmn):
            self.efmn[...] = efmn
        
        @property
        def ofmn(self):
            """
            Element ofmn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__ofmn(self._handle)
            if array_handle in self._arrays:
                ofmn = self._arrays[array_handle]
            else:
                ofmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__ofmn)
                self._arrays[array_handle] = ofmn
            return ofmn
        
        @ofmn.setter
        def ofmn(self, ofmn):
            self.ofmn[...] = ofmn
        
        @property
        def cfmn(self):
            """
            Element cfmn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__cfmn(self._handle)
            if array_handle in self._arrays:
                cfmn = self._arrays[array_handle]
            else:
                cfmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__cfmn)
                self._arrays[array_handle] = cfmn
            return cfmn
        
        @cfmn.setter
        def cfmn(self, cfmn):
            self.cfmn[...] = cfmn
        
        @property
        def sfmn(self):
            """
            Element sfmn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__sfmn(self._handle)
            if array_handle in self._arrays:
                sfmn = self._arrays[array_handle]
            else:
                sfmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__sfmn)
                self._arrays[array_handle] = sfmn
            return sfmn
        
        @sfmn.setter
        def sfmn(self, sfmn):
            self.sfmn[...] = sfmn
        
        @property
        def evmn(self):
            """
            Element evmn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__evmn(self._handle)
            if array_handle in self._arrays:
                evmn = self._arrays[array_handle]
            else:
                evmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__evmn)
                self._arrays[array_handle] = evmn
            return evmn
        
        @evmn.setter
        def evmn(self, evmn):
            self.evmn[...] = evmn
        
        @property
        def odmn(self):
            """
            Element odmn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__odmn(self._handle)
            if array_handle in self._arrays:
                odmn = self._arrays[array_handle]
            else:
                odmn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__odmn)
                self._arrays[array_handle] = odmn
            return odmn
        
        @odmn.setter
        def odmn(self, odmn):
            self.odmn[...] = odmn
        
        @property
        def ijreal(self):
            """
            Element ijreal ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 44
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__ijreal(self._handle)
            if array_handle in self._arrays:
                ijreal = self._arrays[array_handle]
            else:
                ijreal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__ijreal)
                self._arrays[array_handle] = ijreal
            return ijreal
        
        @ijreal.setter
        def ijreal(self, ijreal):
            self.ijreal[...] = ijreal
        
        @property
        def jireal(self):
            """
            Element jireal ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 44
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__jireal(self._handle)
            if array_handle in self._arrays:
                jireal = self._arrays[array_handle]
            else:
                jireal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__jireal)
                self._arrays[array_handle] = jireal
            return jireal
        
        @jireal.setter
        def jireal(self, jireal):
            self.jireal[...] = jireal
        
        @property
        def jkreal(self):
            """
            Element jkreal ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 44
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__jkreal(self._handle)
            if array_handle in self._arrays:
                jkreal = self._arrays[array_handle]
            else:
                jkreal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__jkreal)
                self._arrays[array_handle] = jkreal
            return jkreal
        
        @jkreal.setter
        def jkreal(self, jkreal):
            self.jkreal[...] = jkreal
        
        @property
        def kjreal(self):
            """
            Element kjreal ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 44
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__kjreal(self._handle)
            if array_handle in self._arrays:
                kjreal = self._arrays[array_handle]
            else:
                kjreal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__kjreal)
                self._arrays[array_handle] = kjreal
            return kjreal
        
        @kjreal.setter
        def kjreal(self, kjreal):
            self.kjreal[...] = kjreal
        
        @property
        def bloweremn(self):
            """
            Element bloweremn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 45
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__bloweremn(self._handle)
            if array_handle in self._arrays:
                bloweremn = self._arrays[array_handle]
            else:
                bloweremn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__bloweremn)
                self._arrays[array_handle] = bloweremn
            return bloweremn
        
        @bloweremn.setter
        def bloweremn(self, bloweremn):
            self.bloweremn[...] = bloweremn
        
        @property
        def bloweromn(self):
            """
            Element bloweromn ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 45
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__bloweromn(self._handle)
            if array_handle in self._arrays:
                bloweromn = self._arrays[array_handle]
            else:
                bloweromn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__bloweromn)
                self._arrays[array_handle] = bloweromn
            return bloweromn
        
        @bloweromn.setter
        def bloweromn(self, bloweromn):
            self.bloweromn[...] = bloweromn
        
        @property
        def gbupper(self):
            """
            Element gbupper ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 46
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__gbupper(self._handle)
            if array_handle in self._arrays:
                gbupper = self._arrays[array_handle]
            else:
                gbupper = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__gbupper)
                self._arrays[array_handle] = gbupper
            return gbupper
        
        @gbupper.setter
        def gbupper(self, gbupper):
            self.gbupper[...] = gbupper
        
        @property
        def blower(self):
            """
            Element blower ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 46
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__blower(self._handle)
            if array_handle in self._arrays:
                blower = self._arrays[array_handle]
            else:
                blower = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__blower)
                self._arrays[array_handle] = blower
            return blower
        
        @blower.setter
        def blower(self, blower):
            self.blower[...] = blower
        
        @property
        def basis(self):
            """
            Element basis ftype=real(8) pytype=float
            
            
            Defined at intghs.fpp line 47
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _spec.f90wrap_intghs_workspace__array__basis(self._handle)
            if array_handle in self._arrays:
                basis = self._arrays[array_handle]
            else:
                basis = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _spec.f90wrap_intghs_workspace__array__basis)
                self._arrays[array_handle] = basis
            return basis
        
        @basis.setter
        def basis(self, basis):
            self.basis[...] = basis
        
        def __str__(self):
            ret = ['<intghs_workspace>{\n']
            ret.append('    efmn : ')
            ret.append(repr(self.efmn))
            ret.append(',\n    ofmn : ')
            ret.append(repr(self.ofmn))
            ret.append(',\n    cfmn : ')
            ret.append(repr(self.cfmn))
            ret.append(',\n    sfmn : ')
            ret.append(repr(self.sfmn))
            ret.append(',\n    evmn : ')
            ret.append(repr(self.evmn))
            ret.append(',\n    odmn : ')
            ret.append(repr(self.odmn))
            ret.append(',\n    ijreal : ')
            ret.append(repr(self.ijreal))
            ret.append(',\n    jireal : ')
            ret.append(repr(self.jireal))
            ret.append(',\n    jkreal : ')
            ret.append(repr(self.jkreal))
            ret.append(',\n    kjreal : ')
            ret.append(repr(self.kjreal))
            ret.append(',\n    bloweremn : ')
            ret.append(repr(self.bloweremn))
            ret.append(',\n    bloweromn : ')
            ret.append(repr(self.bloweromn))
            ret.append(',\n    gbupper : ')
            ret.append(repr(self.gbupper))
            ret.append(',\n    blower : ')
            ret.append(repr(self.blower))
            ret.append(',\n    basis : ')
            ret.append(repr(self.basis))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

intghs_module = Intghs_Module()

class Newtontime(f90wrap.runtime.FortranModule):
    """
    Module newtontime
    
    
    Defined at newton.fpp lines 38-40
    
    """
    @property
    def nfcalls(self):
        """
        Element nfcalls ftype=integer  pytype=int
        
        
        Defined at newton.fpp line 39
        
        """
        return _spec.f90wrap_newtontime__get__nfcalls()
    
    @nfcalls.setter
    def nfcalls(self, nfcalls):
        _spec.f90wrap_newtontime__set__nfcalls(nfcalls)
    
    @property
    def ndcalls(self):
        """
        Element ndcalls ftype=integer  pytype=int
        
        
        Defined at newton.fpp line 39
        
        """
        return _spec.f90wrap_newtontime__get__ndcalls()
    
    @ndcalls.setter
    def ndcalls(self, ndcalls):
        _spec.f90wrap_newtontime__set__ndcalls(ndcalls)
    
    @property
    def lastcpu(self):
        """
        Element lastcpu ftype=real(8) pytype=float
        
        
        Defined at newton.fpp line 40
        
        """
        return _spec.f90wrap_newtontime__get__lastcpu()
    
    @lastcpu.setter
    def lastcpu(self, lastcpu):
        _spec.f90wrap_newtontime__set__lastcpu(lastcpu)
    
    def __str__(self):
        ret = ['<newtontime>{\n']
        ret.append('    nfcalls : ')
        ret.append(repr(self.nfcalls))
        ret.append(',\n    ndcalls : ')
        ret.append(repr(self.ndcalls))
        ret.append(',\n    lastcpu : ')
        ret.append(repr(self.lastcpu))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

newtontime = Newtontime()

def preset():
    """
    preset()
    
    
    Defined at preset.fpp lines 16-2667
    
    
    """
    _spec.f90wrap_preset()

def manual():
    """
    manual()
    
    
    Defined at manual.fpp lines 198-230
    
    
    """
    _spec.f90wrap_manual()

def rzaxis(mvol, mn, inrbc, inzbs, inrbs, inzbc, ivol, lcomputederivatives):
    """
    rzaxis(mvol, mn, inrbc, inzbs, inrbs, inzbc, ivol, lcomputederivatives)
    
    
    Defined at rzaxis.fpp lines 60-696
    
    Parameters
    ----------
    mvol : int
    mn : int
    inrbc : float array
    inzbs : float array
    inrbs : float array
    inzbc : float array
    ivol : int
    lcomputederivatives : bool
    
    """
    _spec.f90wrap_rzaxis(mvol=mvol, mn=mn, inrbc=inrbc, inzbs=inzbs, inrbs=inrbs, \
        inzbc=inzbc, ivol=ivol, lcomputederivatives=lcomputederivatives)

def packxi(ngdof, position, mvol, mn, irbc, izbs, irbs, izbc, packorunpack, \
    lcomputederivatives, lcomputeaxis):
    """
    packxi(ngdof, position, mvol, mn, irbc, izbs, irbs, izbc, packorunpack, \
        lcomputederivatives, lcomputeaxis)
    
    
    Defined at packxi.fpp lines 66-165
    
    Parameters
    ----------
    ngdof : int
    position : float array
    mvol : int
    mn : int
    irbc : float array
    izbs : float array
    irbs : float array
    izbc : float array
    packorunpack : str
    lcomputederivatives : bool
    lcomputeaxis : bool
    
    """
    _spec.f90wrap_packxi(ngdof=ngdof, position=position, mvol=mvol, mn=mn, \
        irbc=irbc, izbs=izbs, irbs=irbs, izbc=izbc, packorunpack=packorunpack, \
        lcomputederivatives=lcomputederivatives, lcomputeaxis=lcomputeaxis)

def volume(lvol, vflag):
    """
    volume(lvol, vflag)
    
    
    Defined at volume.fpp lines 31-239
    
    Parameters
    ----------
    lvol : int
    vflag : int
    
    """
    _spec.f90wrap_volume(lvol=lvol, vflag=vflag)

def coords(lvol, lss, lcurvature, ntz, mn):
    """
    coords(lvol, lss, lcurvature, ntz, mn)
    
    
    Defined at coords.fpp lines 115-552
    
    Parameters
    ----------
    lvol : int
    lss : float
    lcurvature : int
    ntz : int
    mn : int
    
    """
    _spec.f90wrap_coords(lvol=lvol, lss=lss, lcurvature=lcurvature, ntz=ntz, mn=mn)

def get_cheby(lss, lrad, cheby):
    """
    get_cheby(lss, lrad, cheby)
    
    
    Defined at basefn.fpp lines 9-35
    
    Parameters
    ----------
    lss : float
    lrad : int
    cheby : float array
    
    """
    _spec.f90wrap_get_cheby(lss=lss, lrad=lrad, cheby=cheby)

def get_cheby_d2(lss, lrad, cheby):
    """
    get_cheby_d2(lss, lrad, cheby)
    
    
    Defined at basefn.fpp lines 37-65
    
    Parameters
    ----------
    lss : float
    lrad : int
    cheby : float array
    
    """
    _spec.f90wrap_get_cheby_d2(lss=lss, lrad=lrad, cheby=cheby)

def get_zernike(r, lrad, mpol, zernike):
    """
    get_zernike(r, lrad, mpol, zernike)
    
    
    Defined at basefn.fpp lines 67-119
    
    Parameters
    ----------
    r : float
    lrad : int
    mpol : int
    zernike : float array
    
    """
    _spec.f90wrap_get_zernike(r=r, lrad=lrad, mpol=mpol, zernike=zernike)

def get_zernike_d2(r, lrad, mpol, zernike):
    """
    get_zernike_d2(r, lrad, mpol, zernike)
    
    
    Defined at basefn.fpp lines 121-181
    
    Parameters
    ----------
    r : float
    lrad : int
    mpol : int
    zernike : float array
    
    """
    _spec.f90wrap_get_zernike_d2(r=r, lrad=lrad, mpol=mpol, zernike=zernike)

def get_zernike_rm(r, lrad, mpol, zernike):
    """
    get_zernike_rm(r, lrad, mpol, zernike)
    
    
    Defined at basefn.fpp lines 183-226
    
    Parameters
    ----------
    r : float
    lrad : int
    mpol : int
    zernike : float array
    
    """
    _spec.f90wrap_get_zernike_rm(r=r, lrad=lrad, mpol=mpol, zernike=zernike)

def allocate_beltrami_matrices(vvol, lcomputederivatives):
    """
    allocate_beltrami_matrices(vvol, lcomputederivatives)
    
    
    Defined at memory.fpp lines 13-126
    
    Parameters
    ----------
    vvol : int
    lcomputederivatives : bool
    
    """
    _spec.f90wrap_allocate_beltrami_matrices(vvol=vvol, \
        lcomputederivatives=lcomputederivatives)

def deallocate_beltrami_matrices(lcomputederivatives):
    """
    deallocate_beltrami_matrices(lcomputederivatives)
    
    
    Defined at memory.fpp lines 129-218
    
    Parameters
    ----------
    lcomputederivatives : bool
    
    """
    \
        _spec.f90wrap_deallocate_beltrami_matrices(lcomputederivatives=lcomputederivatives)

def allocate_geometry_matrices(vvol, lcomputederivatives):
    """
    allocate_geometry_matrices(vvol, lcomputederivatives)
    
    
    Defined at memory.fpp lines 221-587
    
    Parameters
    ----------
    vvol : int
    lcomputederivatives : bool
    
    """
    _spec.f90wrap_allocate_geometry_matrices(vvol=vvol, \
        lcomputederivatives=lcomputederivatives)

def deallocate_geometry_matrices(lcomputederivatives):
    """
    deallocate_geometry_matrices(lcomputederivatives)
    
    
    Defined at memory.fpp lines 590-853
    
    Parameters
    ----------
    lcomputederivatives : bool
    
    """
    \
        _spec.f90wrap_deallocate_geometry_matrices(lcomputederivatives=lcomputederivatives)

def metrix(lquad, lvol):
    """
    metrix(lquad, lvol)
    
    
    Defined at metrix.fpp lines 35-110
    
    Parameters
    ----------
    lquad : int
    lvol : int
    
    """
    _spec.f90wrap_metrix(lquad=lquad, lvol=lvol)

def compute_guvijsave(lquad, vvol, ideriv, lcurvature):
    """
    compute_guvijsave(lquad, vvol, ideriv, lcurvature)
    
    
    Defined at metrix.fpp lines 113-129
    
    Parameters
    ----------
    lquad : int
    vvol : int
    ideriv : int
    lcurvature : int
    
    """
    _spec.f90wrap_compute_guvijsave(lquad=lquad, vvol=vvol, ideriv=ideriv, \
        lcurvature=lcurvature)

def ma00aa(lquad, mn, lvol, lrad):
    """
    ma00aa(lquad, mn, lvol, lrad)
    
    
    Defined at ma00aa.fpp lines 40-317
    
    Parameters
    ----------
    lquad : int
    mn : int
    lvol : int
    lrad : int
    
    """
    _spec.f90wrap_ma00aa(lquad=lquad, mn=mn, lvol=lvol, lrad=lrad)

def matrix(lvol, mn, lrad):
    """
    matrix(lvol, mn, lrad)
    
    
    Defined at matrix.fpp lines 232-499
    
    Parameters
    ----------
    lvol : int
    mn : int
    lrad : int
    
    """
    _spec.f90wrap_matrix(lvol=lvol, mn=mn, lrad=lrad)

def matrixbg(lvol, mn, lrad):
    """
    matrixbg(lvol, mn, lrad)
    
    
    Defined at matrix.fpp lines 502-533
    
    Parameters
    ----------
    lvol : int
    mn : int
    lrad : int
    
    """
    _spec.f90wrap_matrixbg(lvol=lvol, mn=mn, lrad=lrad)

def spsmat(lvol, mn, lrad):
    """
    spsmat(lvol, mn, lrad)
    
    
    Defined at spsmat.fpp lines 12-421
    
    Parameters
    ----------
    lvol : int
    mn : int
    lrad : int
    
    """
    _spec.f90wrap_spsmat(lvol=lvol, mn=mn, lrad=lrad)

def push_back(iq, nq, nn, va, vd, vja, qa, qd, qja):
    """
    push_back(iq, nq, nn, va, vd, vja, qa, qd, qja)
    
    
    Defined at spsmat.fpp lines 425-446
    
    Parameters
    ----------
    iq : int
    nq : int array
    nn : int
    va : float
    vd : float
    vja : int
    qa : float array
    qd : float array
    qja : int array
    
    """
    _spec.f90wrap_push_back(iq=iq, nq=nq, nn=nn, va=va, vd=vd, vja=vja, qa=qa, \
        qd=qd, qja=qja)

def clean_queue(nq, nn, qa, qd, qja):
    """
    clean_queue(nq, nn, qa, qd, qja)
    
    
    Defined at spsmat.fpp lines 448-459
    
    Parameters
    ----------
    nq : int array
    nn : int
    qa : float array
    qd : float array
    qja : int array
    
    """
    _spec.f90wrap_clean_queue(nq=nq, nn=nn, qa=qa, qd=qd, qja=qja)

def addline(nq, nn, qa, qd, qja, ns, nrow, dmas, dmds, jdmas, idmas):
    """
    addline(nq, nn, qa, qd, qja, ns, nrow, dmas, dmds, jdmas, idmas)
    
    
    Defined at spsmat.fpp lines 461-477
    
    Parameters
    ----------
    nq : int array
    nn : int
    qa : float array
    qd : float array
    qja : int array
    ns : int
    nrow : int
    dmas : float array
    dmds : float array
    jdmas : int array
    idmas : int array
    
    """
    _spec.f90wrap_addline(nq=nq, nn=nn, qa=qa, qd=qd, qja=qja, ns=ns, nrow=nrow, \
        dmas=dmas, dmds=dmds, jdmas=jdmas, idmas=idmas)

def spsint(lquad, mn, lvol, lrad):
    """
    spsint(lquad, mn, lvol, lrad)
    
    
    Defined at spsint.fpp lines 12-220
    
    Parameters
    ----------
    lquad : int
    mn : int
    lvol : int
    lrad : int
    
    """
    _spec.f90wrap_spsint(lquad=lquad, mn=mn, lvol=lvol, lrad=lrad)

def mp00ac(ndof, xdof, fdof, ddof, ldfjac, iflag):
    """
    mp00ac(ndof, xdof, fdof, ddof, ldfjac, iflag)
    
    
    Defined at mp00ac.fpp lines 89-811
    
    Parameters
    ----------
    ndof : int
    xdof : float array
    fdof : float array
    ddof : float array
    ldfjac : int
    iflag : int
    
    """
    _spec.f90wrap_mp00ac(ndof=ndof, xdof=xdof, fdof=fdof, ddof=ddof, ldfjac=ldfjac, \
        iflag=iflag)

def rungmres(n, nrestart, mu, vvol, rhs, sol, ipar, fpar, wk, nw, guess, a, au, \
    jau, ju, iperm, ierr):
    """
    rungmres(n, nrestart, mu, vvol, rhs, sol, ipar, fpar, wk, nw, guess, a, au, jau, \
        ju, iperm, ierr)
    
    
    Defined at mp00ac.fpp lines 814-868
    
    Parameters
    ----------
    n : int
    nrestart : int
    mu : float
    vvol : int
    rhs : float array
    sol : float array
    ipar : int array
    fpar : float array
    wk : float array
    nw : int
    guess : float array
    a : float array
    au : float array
    jau : int array
    ju : int array
    iperm : int array
    ierr : int
    
    """
    _spec.f90wrap_rungmres(n=n, nrestart=nrestart, mu=mu, vvol=vvol, rhs=rhs, \
        sol=sol, ipar=ipar, fpar=fpar, wk=wk, nw=nw, guess=guess, a=a, au=au, \
        jau=jau, ju=ju, iperm=iperm, ierr=ierr)

def matvec(n, x, ax, a, mu, vvol):
    """
    matvec(n, x, ax, a, mu, vvol)
    
    
    Defined at mp00ac.fpp lines 870-893
    
    Parameters
    ----------
    n : int
    x : float array
    ax : float array
    a : float array
    mu : float
    vvol : int
    
    """
    _spec.f90wrap_matvec(n=n, x=x, ax=ax, a=a, mu=mu, vvol=vvol)

def prec_solve(n, vecin, vecout, au, jau, ju, iperm):
    """
    prec_solve(n, vecin, vecout, au, jau, ju, iperm)
    
    
    Defined at mp00ac.fpp lines 895-907
    
    Parameters
    ----------
    n : int
    vecin : float array
    vecout : float array
    au : float array
    jau : int array
    ju : int array
    iperm : int array
    
    """
    _spec.f90wrap_prec_solve(n=n, vecin=vecin, vecout=vecout, au=au, jau=jau, ju=ju, \
        iperm=iperm)

def ma02aa(lvol, nn):
    """
    ma02aa(lvol, nn)
    
    
    Defined at ma02aa.fpp lines 16-673
    
    Parameters
    ----------
    lvol : int
    nn : int
    
    """
    _spec.f90wrap_ma02aa(lvol=lvol, nn=nn)

def packab(packorunpack, lvol, nn, solution, ideriv):
    """
    packab(packorunpack, lvol, nn, solution, ideriv)
    
    
    Defined at packab.fpp lines 28-197
    
    Parameters
    ----------
    packorunpack : str
    lvol : int
    nn : int
    solution : float array
    ideriv : int
    
    """
    _spec.f90wrap_packab(packorunpack=packorunpack, lvol=lvol, nn=nn, \
        solution=solution, ideriv=ideriv)

def tr00ab(lvol, mn, nn, nt, nz, iflag, ldiota):
    """
    tr00ab(lvol, mn, nn, nt, nz, iflag, ldiota)
    
    
    Defined at tr00ab.fpp lines 58-503
    
    Parameters
    ----------
    lvol : int
    mn : int
    nn : int
    nt : int
    nz : int
    iflag : int
    ldiota : float array
    
    """
    _spec.f90wrap_tr00ab(lvol=lvol, mn=mn, nn=nn, nt=nt, nz=nz, iflag=iflag, \
        ldiota=ldiota)

def curent(lvol, mn, nt, nz, iflag, lditgp):
    """
    curent(lvol, mn, nt, nz, iflag, lditgp)
    
    
    Defined at curent.fpp lines 54-163
    
    Parameters
    ----------
    lvol : int
    mn : int
    nt : int
    nz : int
    iflag : int
    lditgp : float array
    
    """
    _spec.f90wrap_curent(lvol=lvol, mn=mn, nt=nt, nz=nz, iflag=iflag, lditgp=lditgp)

def df00ab(pnn, xi, fxi, dfxi, ldfjac, iflag):
    """
    df00ab(pnn, xi, fxi, dfxi, ldfjac, iflag)
    
    
    Defined at df00ab.fpp lines 16-82
    
    Parameters
    ----------
    pnn : int
    xi : float array
    fxi : float array
    dfxi : float array
    ldfjac : int
    iflag : int
    
    """
    _spec.f90wrap_df00ab(pnn=pnn, xi=xi, fxi=fxi, dfxi=dfxi, ldfjac=ldfjac, \
        iflag=iflag)

def lforce(lvol, iocons, ideriv, ntz, dbb, xx, yy, length, ddl, mml, iflag):
    """
    lforce(lvol, iocons, ideriv, ntz, dbb, xx, yy, length, ddl, mml, iflag)
    
    
    Defined at lforce.fpp lines 125-295
    
    Parameters
    ----------
    lvol : int
    iocons : int
    ideriv : int
    ntz : int
    dbb : float array
    xx : float array
    yy : float array
    length : float array
    ddl : float
    mml : float
    iflag : int
    
    """
    _spec.f90wrap_lforce(lvol=lvol, iocons=iocons, ideriv=ideriv, ntz=ntz, dbb=dbb, \
        xx=xx, yy=yy, length=length, ddl=ddl, mml=mml, iflag=iflag)

def intghs(lquad, mn, lvol, lrad, idx):
    """
    intghs(lquad, mn, lvol, lrad, idx)
    
    
    Defined at intghs.fpp lines 52-252
    
    Parameters
    ----------
    lquad : int
    mn : int
    lvol : int
    lrad : int
    idx : int
    
    """
    _spec.f90wrap_intghs(lquad=lquad, mn=mn, lvol=lvol, lrad=lrad, idx=idx)

def intghs_workspace_init(lvol):
    """
    intghs_workspace_init(lvol)
    
    
    Defined at intghs.fpp lines 255-404
    
    Parameters
    ----------
    lvol : int
    
    """
    _spec.f90wrap_intghs_workspace_init(lvol=lvol)

def intghs_workspace_destroy():
    """
    intghs_workspace_destroy()
    
    
    Defined at intghs.fpp lines 406-520
    
    
    """
    _spec.f90wrap_intghs_workspace_destroy()

def mtrxhs(lvol, mn, lrad, resulta, resultd, idx):
    """
    mtrxhs(lvol, mn, lrad, resulta, resultd, idx)
    
    
    Defined at mtrxhs.fpp lines 12-204
    
    Parameters
    ----------
    lvol : int
    mn : int
    lrad : int
    resulta : float array
    resultd : float array
    idx : int
    
    """
    _spec.f90wrap_mtrxhs(lvol=lvol, mn=mn, lrad=lrad, resulta=resulta, \
        resultd=resultd, idx=idx)

def lbpol(lvol, bt00, ideriv, iocons):
    """
    lbpol(lvol, bt00, ideriv, iocons)
    
    
    Defined at lbpol.fpp lines 28-127
    
    Parameters
    ----------
    lvol : int
    bt00 : float array
    ideriv : int
    iocons : int
    
    """
    _spec.f90wrap_lbpol(lvol=lvol, bt00=bt00, ideriv=ideriv, iocons=iocons)

def brcast(lvol):
    """
    brcast(lvol)
    
    
    Defined at brcast.fpp lines 25-228
    
    Parameters
    ----------
    lvol : int
    
    """
    _spec.f90wrap_brcast(lvol=lvol)

def dfp100(ndofgl, x, fvec, lcomputederivatives):
    """
    dfp100(ndofgl, x, fvec, lcomputederivatives)
    
    
    Defined at dfp100.fpp lines 31-311
    
    Parameters
    ----------
    ndofgl : int
    x : float array
    fvec : float array
    lcomputederivatives : bool
    
    ------
     vvol:                       loop index on volumes
     Ndofgl: Input parameter necessary for the use of hybrd1. Unused otherwise.
     iflag:                      Flag changed by hybrd1
     cpu_send_one, cpu_send_two: CPU IDs, used for MPI communications
     status:                     MPI status
     Fvec:                       Global constraint values
     x: Degrees of freedom of hybrd1. For now contains only the poloidal flux
    """
    _spec.f90wrap_dfp100(ndofgl=ndofgl, x=x, fvec=fvec, \
        lcomputederivatives=lcomputederivatives)

def dfp200(lcomputederivatives, vvol):
    """
    dfp200(lcomputederivatives, vvol)
    
    
    Defined at dfp200.fpp lines 38-781
    
    Parameters
    ----------
    lcomputederivatives : bool
    vvol : int
    
    """
    _spec.f90wrap_dfp200(lcomputederivatives=lcomputederivatives, vvol=vvol)

def get_lu_beltrami_matrices(vvol, obi, nn):
    """
    get_lu_beltrami_matrices(vvol, obi, nn)
    
    
    Defined at dfp200.fpp lines 786-891
    
    Parameters
    ----------
    vvol : int
    obi : Matrixlu
    nn : int
    
    """
    _spec.f90wrap_get_lu_beltrami_matrices(vvol=vvol, obi=obi._handle, nn=nn)

def get_perturbed_solution(vvol, obi, nn):
    """
    get_perturbed_solution(vvol, obi, nn)
    
    
    Defined at dfp200.fpp lines 894-961
    
    Parameters
    ----------
    vvol : int
    obi : Matrixlu
    nn : int
    
    ------
    """
    _spec.f90wrap_get_perturbed_solution(vvol=vvol, obi=obi._handle, nn=nn)

def evaluate_dmupfdx(innout, idof, ii, issym, irz):
    """
    evaluate_dmupfdx(innout, idof, ii, issym, irz)
    
    
    Defined at dfp200.fpp lines 964-1303
    
    Parameters
    ----------
    innout : int
    idof : int
    ii : int
    issym : int
    irz : int
    
    """
    _spec.f90wrap_evaluate_dmupfdx(innout=innout, idof=idof, ii=ii, issym=issym, \
        irz=irz)

def evaluate_dbb(lvol, idof, innout, issym, irz, ii, dbb, xx, yy, length, drr, \
    dzz, dii, dll, dpp, ntz):
    """
    evaluate_dbb(lvol, idof, innout, issym, irz, ii, dbb, xx, yy, length, drr, dzz, \
        dii, dll, dpp, ntz)
    
    
    Defined at dfp200.fpp lines 1306-1615
    
    Parameters
    ----------
    lvol : int
    idof : int
    innout : int
    issym : int
    irz : int
    ii : int
    dbb : float array
    xx : float array
    yy : float array
    length : float array
    drr : float array
    dzz : float array
    dii : float array
    dll : float array
    dpp : float array
    ntz : int
    
    ------
    """
    _spec.f90wrap_evaluate_dbb(lvol=lvol, idof=idof, innout=innout, issym=issym, \
        irz=irz, ii=ii, dbb=dbb, xx=xx, yy=yy, length=length, drr=drr, dzz=dzz, \
        dii=dii, dll=dll, dpp=dpp, ntz=ntz)

def dforce(ngdof, position, force, lcomputederivatives, lcomputeaxis):
    """
    dforce(ngdof, position, force, lcomputederivatives, lcomputeaxis)
    
    
    Defined at dforce.fpp lines 80-656
    
    Parameters
    ----------
    ngdof : int
    position : float array
    force : float array
    lcomputederivatives : bool
    lcomputeaxis : bool
    
    """
    _spec.f90wrap_dforce(ngdof=ngdof, position=position, force=force, \
        lcomputederivatives=lcomputederivatives, lcomputeaxis=lcomputeaxis)

def newton(ngdof, position):
    """
    ihybrd = newton(ngdof, position)
    
    
    Defined at newton.fpp lines 43-356
    
    Parameters
    ----------
    ngdof : int
    position : float array
    
    Returns
    -------
    ihybrd : int
    
    """
    ihybrd = _spec.f90wrap_newton(ngdof=ngdof, position=position)
    return ihybrd

def writereadgf(readorwrite, ngdof):
    """
    ireadhessian = writereadgf(readorwrite, ngdof)
    
    
    Defined at newton.fpp lines 360-485
    
    Parameters
    ----------
    readorwrite : str
    ngdof : int
    
    Returns
    -------
    ireadhessian : int
    
    """
    ireadhessian = _spec.f90wrap_writereadgf(readorwrite=readorwrite, ngdof=ngdof)
    return ireadhessian

def fcn1(ngdof, xx, fvec, irevcm):
    """
    fcn1(ngdof, xx, fvec, irevcm)
    
    
    Defined at newton.fpp lines 489-615
    
    Parameters
    ----------
    ngdof : int
    xx : float array
    fvec : float array
    irevcm : int
    
    """
    _spec.f90wrap_fcn1(ngdof=ngdof, xx=xx, fvec=fvec, irevcm=irevcm)

def fcn2(ngdof, xx, fvec, fjac, ldfjac, irevcm):
    """
    fcn2(ngdof, xx, fvec, fjac, ldfjac, irevcm)
    
    
    Defined at newton.fpp lines 619-782
    
    Parameters
    ----------
    ngdof : int
    xx : float array
    fvec : float array
    fjac : float array
    ldfjac : int
    irevcm : int
    
    """
    _spec.f90wrap_fcn2(ngdof=ngdof, xx=xx, fvec=fvec, fjac=fjac, ldfjac=ldfjac, \
        irevcm=irevcm)

def casing(teta, zeta, icasing):
    """
    gbn = casing(teta, zeta, icasing)
    
    
    Defined at casing.fpp lines 77-207
    
    Parameters
    ----------
    teta : float
    zeta : float
    icasing : int
    
    Returns
    -------
    gbn : float
    
    """
    gbn = _spec.f90wrap_casing(teta=teta, zeta=zeta, icasing=icasing)
    return gbn

def dvcfield(ndim, tz, nfun, vcintegrand):
    """
    dvcfield(ndim, tz, nfun, vcintegrand)
    
    
    Defined at casing.fpp lines 234-442
    
    Parameters
    ----------
    ndim : int
    tz : float array
    nfun : int
    vcintegrand : float array
    
    """
    _spec.f90wrap_dvcfield(ndim=ndim, tz=tz, nfun=nfun, vcintegrand=vcintegrand)

def bnorml(mn, ntz, efmn, ofmn):
    """
    bnorml(mn, ntz, efmn, ofmn)
    
    
    Defined at bnorml.fpp lines 61-289
    
    Parameters
    ----------
    mn : int
    ntz : int
    efmn : float array
    ofmn : float array
    
    """
    _spec.f90wrap_bnorml(mn=mn, ntz=ntz, efmn=efmn, ofmn=ofmn)

def vcintegrand(lteta, lzeta):
    """
    vcintegrand = vcintegrand(lteta, lzeta)
    
    
    Defined at bnorml.fpp lines 370-559
    
    Parameters
    ----------
    lteta : float
    lzeta : float
    
    Returns
    -------
    vcintegrand : float
    
    """
    vcintegrand = _spec.f90wrap_vcintegrand(lteta=lteta, lzeta=lzeta)
    return vcintegrand

def zetalow(teta):
    """
    zetalow = zetalow(teta)
    
    
    Defined at bnorml.fpp lines 563-574
    
    Parameters
    ----------
    teta : float
    
    Returns
    -------
    zetalow : float
    
    """
    zetalow = _spec.f90wrap_zetalow(teta=teta)
    return zetalow

def zetaupp(teta):
    """
    zetaupp = zetaupp(teta)
    
    
    Defined at bnorml.fpp lines 578-589
    
    Parameters
    ----------
    teta : float
    
    Returns
    -------
    zetaupp : float
    
    """
    zetaupp = _spec.f90wrap_zetaupp(teta=teta)
    return zetaupp

def jo00aa(lvol, ntz, lquad, mn):
    """
    jo00aa(lvol, ntz, lquad, mn)
    
    
    Defined at jo00aa.fpp lines 46-374
    
    Parameters
    ----------
    lvol : int
    ntz : int
    lquad : int
    mn : int
    
    """
    _spec.f90wrap_jo00aa(lvol=lvol, ntz=ntz, lquad=lquad, mn=mn)

def pp00aa():
    """
    pp00aa()
    
    
    Defined at pp00aa.fpp lines 69-306
    
    
    """
    _spec.f90wrap_pp00aa()

def pp00ab(lvol, sti, nz, nppts, poincaredata, fittedtransform):
    """
    utflag = pp00ab(lvol, sti, nz, nppts, poincaredata, fittedtransform)
    
    
    Defined at pp00ab.fpp lines 33-159
    
    Parameters
    ----------
    lvol : int
    sti : float array
    nz : int
    nppts : int
    poincaredata : float array
    fittedtransform : float array
    
    Returns
    -------
    utflag : int
    
    """
    utflag = _spec.f90wrap_pp00ab(lvol=lvol, sti=sti, nz=nz, nppts=nppts, \
        poincaredata=poincaredata, fittedtransform=fittedtransform)
    return utflag

def bfield(zeta, st, bst):
    """
    bfield(zeta, st, bst)
    
    
    Defined at bfield.fpp lines 30-137
    
    Parameters
    ----------
    zeta : float
    st : float array
    bst : float array
    
    """
    _spec.f90wrap_bfield(zeta=zeta, st=st, bst=bst)

def bfield_tangent(zeta, st, bst):
    """
    bfield_tangent(zeta, st, bst)
    
    
    Defined at bfield.fpp lines 140-277
    
    Parameters
    ----------
    zeta : float
    st : float array
    bst : float array
    
    """
    _spec.f90wrap_bfield_tangent(zeta=zeta, st=st, bst=bst)

def stzxyz(lvol, stz, rpz):
    """
    stzxyz(lvol, stz, rpz)
    
    
    Defined at stzxyz.fpp lines 22-119
    
    Parameters
    ----------
    lvol : int
    stz : float array
    rpz : float array
    
    """
    _spec.f90wrap_stzxyz(lvol=lvol, stz=stz, rpz=rpz)

def hesian(ngdof, position, mvol, mn, lgdof):
    """
    hesian(ngdof, position, mvol, mn, lgdof)
    
    
    Defined at hesian.fpp lines 17-556
    
    Parameters
    ----------
    ngdof : int
    position : float array
    mvol : int
    mn : int
    lgdof : int
    
    """
    _spec.f90wrap_hesian(ngdof=ngdof, position=position, mvol=mvol, mn=mn, \
        lgdof=lgdof)

def ra00aa(writeorread):
    """
    ra00aa(writeorread)
    
    
    Defined at ra00aa.fpp lines 37-282
    
    Parameters
    ----------
    writeorread : str
    
    """
    _spec.f90wrap_ra00aa(writeorread=writeorread)

def gi00ab(mpol, ntor, nfp, mn, im, in_):
    """
    gi00ab(mpol, ntor, nfp, mn, im, in_)
    
    
    Defined at numrec.fpp lines 48-64
    
    Parameters
    ----------
    mpol : int
    ntor : int
    nfp : int
    mn : int
    im : int array
    in_ : int array
    
    """
    _spec.f90wrap_gi00ab(mpol=mpol, ntor=ntor, nfp=nfp, mn=mn, im=im, in_=in_)

def getimn(mpol, ntor, nfp, mi, ni):
    """
    idx = getimn(mpol, ntor, nfp, mi, ni)
    
    
    Defined at numrec.fpp lines 67-78
    
    Parameters
    ----------
    mpol : int
    ntor : int
    nfp : int
    mi : int
    ni : int
    
    Returns
    -------
    idx : int
    
    """
    idx = _spec.f90wrap_getimn(mpol=mpol, ntor=ntor, nfp=nfp, mi=mi, ni=ni)
    return idx

def tfft(nt, nz, ijreal, ijimag, mn, im, in_, efmn, ofmn, cfmn, sfmn, ifail):
    """
    tfft(nt, nz, ijreal, ijimag, mn, im, in_, efmn, ofmn, cfmn, sfmn, ifail)
    
    
    Defined at numrec.fpp lines 94-140
    
    Parameters
    ----------
    nt : int
    nz : int
    ijreal : float array
    ijimag : float array
    mn : int
    im : int array
    in_ : int array
    efmn : float array
    ofmn : float array
    cfmn : float array
    sfmn : float array
    ifail : int
    
    """
    _spec.f90wrap_tfft(nt=nt, nz=nz, ijreal=ijreal, ijimag=ijimag, mn=mn, im=im, \
        in_=in_, efmn=efmn, ofmn=ofmn, cfmn=cfmn, sfmn=sfmn, ifail=ifail)

def invfft(mn, im, in_, efmn, ofmn, cfmn, sfmn, nt, nz, ijreal, ijimag):
    """
    invfft(mn, im, in_, efmn, ofmn, cfmn, sfmn, nt, nz, ijreal, ijimag)
    
    
    Defined at numrec.fpp lines 148-176
    
    Parameters
    ----------
    mn : int
    im : int array
    in_ : int array
    efmn : float array
    ofmn : float array
    cfmn : float array
    sfmn : float array
    nt : int
    nz : int
    ijreal : float array
    ijimag : float array
    
    """
    _spec.f90wrap_invfft(mn=mn, im=im, in_=in_, efmn=efmn, ofmn=ofmn, cfmn=cfmn, \
        sfmn=sfmn, nt=nt, nz=nz, ijreal=ijreal, ijimag=ijimag)

def gauleg(n, weight, abscis):
    """
    ifail = gauleg(n, weight, abscis)
    
    
    Defined at numrec.fpp lines 184-224
    
    Parameters
    ----------
    n : int
    weight : float array
    abscis : float array
    
    Returns
    -------
    ifail : int
    
    """
    ifail = _spec.f90wrap_gauleg(n=n, weight=weight, abscis=abscis)
    return ifail

def read_command_args():
    """
    read_command_args()
    
    
    Defined at xspech.fpp lines 150-206
    
    
    """
    _spec.f90wrap_read_command_args()

def spec():
    """
    spec()
    
    
    Defined at xspech.fpp lines 210-671
    
    
    """
    _spec.f90wrap_spec()

def final_diagnostics():
    """
    final_diagnostics()
    
    
    Defined at xspech.fpp lines 681-856
    
    
    """
    _spec.f90wrap_final_diagnostics()

def ending():
    """
    ending()
    
    
    Defined at xspech.fpp lines 860-1196
    
    
    """
    _spec.f90wrap_ending()

