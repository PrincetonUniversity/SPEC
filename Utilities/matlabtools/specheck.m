function matching = specheck(fdata, gdata, idata, pdata, data)

% Compares the outputs of the current SPEC version and the newly-developed single HDF5 output.
%
% INPUT
% -  fdata    : from read_spec_field('master/testcase.sp.h5')
% -  gdata    : from read_spec_grid('master/testcase.sp.h5')
% -  idata    : from read_spec_iota('master/testcase.sp.h5')
% -  pdata    : from read_spec_poincare('master/testcase.sp.h5')
% -  data     : from read_spec('issue68/testcase.h5')
%
% OUTPUT
% -  matching : 1 if the outputs are identical, 0 otherwise
%
% written by J.Schilling (2019)

matching = 1;

% check equivalence of fdata: from "testcase.sp.h5" file
% quantities from old .sp.h5 file, sorted by equivalent quantity in the new output file
if isequal(fdata.forcetol    , data.input.global.forcetol       ) ; disp('ok: fdata.forcetol   ') ; else disp('ERROR: fdata.forcetol   '); matching=0; end

if isequal(fdata.Igeometry   , data.input.physics.Igeometry     ) ; disp('ok: fdata.Igeometry  ') ; else disp('ERROR: fdata.Igeometry  '); matching=0; end
if isequal(fdata.Istellsym   , data.input.physics.Istellsym     ) ; disp('ok: fdata.Istellsym  ') ; else disp('ERROR: fdata.Istellsym  '); matching=0; end
if isequal(fdata.Ladiabatic  , data.input.physics.Ladiabatic    ) ; disp('ok: fdata.Ladiabatic ') ; else disp('ERROR: fdata.Ladiabatic '); matching=0; end
if isequal(fdata.Lconstraint , data.input.physics.Lconstraint   ) ; disp('ok: fdata.Lconstraint') ; else disp('ERROR: fdata.Lconstraint'); matching=0; end
if isequal(fdata.Lfreebound  , data.input.physics.Lfreebound    ) ; disp('ok: fdata.Lfreebound ') ; else disp('ERROR: fdata.Lfreebound '); matching=0; end
if isequal(fdata.Lrad        , data.input.physics.Lrad          ) ; disp('ok: fdata.Lrad       ') ; else disp('ERROR: fdata.Lrad       '); matching=0; end
if isequal(fdata.Mpol        , data.input.physics.Mpol          ) ; disp('ok: fdata.Mpol       ') ; else disp('ERROR: fdata.Mpol       '); matching=0; end
if isequal(fdata.Nfp         , data.input.physics.Nfp           ) ; disp('ok: fdata.Nfp        ') ; else disp('ERROR: fdata.Nfp        '); matching=0; end
if isequal(fdata.Ntor        , data.input.physics.Ntor          ) ; disp('ok: fdata.Ntor       ') ; else disp('ERROR: fdata.Ntor       '); matching=0; end
if isequal(fdata.Nvol        , data.input.physics.Nvol          ) ; disp('ok: fdata.Nvol       ') ; else disp('ERROR: fdata.Nvol       '); matching=0; end
if isequal(fdata.curpol      , data.input.physics.curpol        ) ; disp('ok: fdata.curpol     ') ; else disp('ERROR: fdata.curpol     '); matching=0; end
if isequal(fdata.curtor      , data.input.physics.curtor        ) ; disp('ok: fdata.curtor     ') ; else disp('ERROR: fdata.curtor     '); matching=0; end
if isequal(fdata.gamma       , data.input.physics.gamma         ) ; disp('ok: fdata.gamma      ') ; else disp('ERROR: fdata.gamma      '); matching=0; end
if isequal(fdata.iota        , data.input.physics.iota          ) ; disp('ok: fdata.iota       ') ; else disp('ERROR: fdata.iota       '); matching=0; end
if isequal(fdata.lp          , data.input.physics.lp            ) ; disp('ok: fdata.lp         ') ; else disp('ERROR: fdata.lp         '); matching=0; end
if isequal(fdata.lq          , data.input.physics.lq            ) ; disp('ok: fdata.lq         ') ; else disp('ERROR: fdata.lq         '); matching=0; end
if isequal(fdata.mu          , data.input.physics.mu            ) ; disp('ok: fdata.mu         ') ; else disp('ERROR: fdata.mu         '); matching=0; end
if isequal(fdata.mupfits     , data.input.physics.mupfits       ) ; disp('ok: fdata.mupfits    ') ; else disp('ERROR: fdata.mupfits    '); matching=0; end
if isequal(fdata.mupftol     , data.input.physics.mupftol       ) ; disp('ok: fdata.mupftol    ') ; else disp('ERROR: fdata.mupftol    '); matching=0; end
if isequal(fdata.oita        , data.input.physics.oita          ) ; disp('ok: fdata.oita       ') ; else disp('ERROR: fdata.oita       '); matching=0; end
if isequal(fdata.pflux       , data.input.physics.pflux         ) ; disp('ok: fdata.pflux      ') ; else disp('ERROR: fdata.pflux      '); matching=0; end
if isequal(fdata.phiedge     , data.input.physics.phiedge       ) ; disp('ok: fdata.phiedge    ') ; else disp('ERROR: fdata.phiedge    '); matching=0; end
if isequal(fdata.pl          , data.input.physics.pl            ) ; disp('ok: fdata.pl         ') ; else disp('ERROR: fdata.pl         '); matching=0; end
if isequal(fdata.pr          , data.input.physics.pr            ) ; disp('ok: fdata.pr         ') ; else disp('ERROR: fdata.pr         '); matching=0; end
if isequal(fdata.pressure    , data.input.physics.pressure      ) ; disp('ok: fdata.pressure   ') ; else disp('ERROR: fdata.pressure   '); matching=0; end
if isequal(fdata.pscale      , data.input.physics.pscale        ) ; disp('ok: fdata.pscale     ') ; else disp('ERROR: fdata.pscale     '); matching=0; end
if isequal(fdata.ql          , data.input.physics.ql            ) ; disp('ok: fdata.ql         ') ; else disp('ERROR: fdata.ql         '); matching=0; end
if isequal(fdata.qr          , data.input.physics.qr            ) ; disp('ok: fdata.qr         ') ; else disp('ERROR: fdata.qr         '); matching=0; end
if isequal(fdata.rp          , data.input.physics.rp            ) ; disp('ok: fdata.rp         ') ; else disp('ERROR: fdata.rp         '); matching=0; end
if isequal(fdata.rq          , data.input.physics.rq            ) ; disp('ok: fdata.rq         ') ; else disp('ERROR: fdata.rq         '); matching=0; end
if isequal(fdata.tflux       , data.input.physics.tflux         ) ; disp('ok: fdata.tflux      ') ; else disp('ERROR: fdata.tflux      '); matching=0; end

if isequal(fdata.Lperturbed  , data.input.diagnostics.Lperturbed) ; disp('ok: fdata.Lperturbed ') ; else disp('ERROR: fdata.Lperturbed '); matching=0; end
if isequal(fdata.dpp         , data.input.diagnostics.dpp       ) ; disp('ok: fdata.dpp        ') ; else disp('ERROR: fdata.dpp        '); matching=0; end
if isequal(fdata.dqq         , data.input.diagnostics.dqq       ) ; disp('ok: fdata.dqq        ') ; else disp('ERROR: fdata.dqq        '); matching=0; end

if isequal(fdata.Bnc         , data.output.Bnc                  ) ; disp('ok: fdata.Bnc        ') ; else disp('ERROR: fdata.Bnc        '); matching=0; end
if isequal(fdata.Bns         , data.output.Bns                  ) ; disp('ok: fdata.Bns        ') ; else disp('ERROR: fdata.Bns        '); matching=0; end
if isequal(fdata.Btemn       , data.output.Btemn                ) ; disp('ok: fdata.Btemn      ') ; else disp('ERROR: fdata.Btemn      '); matching=0; end
if isequal(fdata.Btomn       , data.output.Btomn                ) ; disp('ok: fdata.Btomn      ') ; else disp('ERROR: fdata.Btomn      '); matching=0; end
if isequal(fdata.Bzemn       , data.output.Bzemn                ) ; disp('ok: fdata.Bzemn      ') ; else disp('ERROR: fdata.Bzemn      '); matching=0; end
if isequal(fdata.Bzomn       , data.output.Bzomn                ) ; disp('ok: fdata.Bzomn      ') ; else disp('ERROR: fdata.Bzomn      '); matching=0; end
if isequal(fdata.ForceErr    , data.output.ForceErr             ) ; disp('ok: fdata.ForceErr   ') ; else disp('ERROR: fdata.ForceErr   '); matching=0; end
if isequal(fdata.Mrad        , data.output.Mrad                 ) ; disp('ok: fdata.Mrad       ') ; else disp('ERROR: fdata.Mrad       '); matching=0; end
if isequal(fdata.Mvol        , data.output.Mvol                 ) ; disp('ok: fdata.Mvol       ') ; else disp('ERROR: fdata.Mvol       '); matching=0; end
if isequal(fdata.Rbc         , data.output.Rbc                  ) ; disp('ok: fdata.Rbc        ') ; else disp('ERROR: fdata.Rbc        '); matching=0; end
if isequal(fdata.Rbs         , data.output.Rbs                  ) ; disp('ok: fdata.Rbs        ') ; else disp('ERROR: fdata.Rbs        '); matching=0; end
if isequal(fdata.TT          , data.output.TT                   ) ; disp('ok: fdata.TT         ') ; else disp('ERROR: fdata.TT         '); matching=0; end
if isequal(fdata.Vnc         , data.output.Vnc                  ) ; disp('ok: fdata.Vnc        ') ; else disp('ERROR: fdata.Vnc        '); matching=0; end
if isequal(fdata.Vns         , data.output.Vns                  ) ; disp('ok: fdata.Vns        ') ; else disp('ERROR: fdata.Vns        '); matching=0; end
if isequal(fdata.Zbc         , data.output.Zbc                  ) ; disp('ok: fdata.Zbc        ') ; else disp('ERROR: fdata.Zbc        '); matching=0; end
if isequal(fdata.Zbs         , data.output.Zbs                  ) ; disp('ok: fdata.Zbs        ') ; else disp('ERROR: fdata.Zbs        '); matching=0; end
if isequal(fdata.adiabatic   , data.output.adiabatic            ) ; disp('ok: fdata.adiabatic  ') ; else disp('ERROR: fdata.adiabatic  '); matching=0; end
if isequal(fdata.helicity    , data.output.helicity             ) ; disp('ok: fdata.helicity   ') ; else disp('ERROR: fdata.helicity   '); matching=0; end
if isequal(fdata.im          , data.output.im                   ) ; disp('ok: fdata.im         ') ; else disp('ERROR: fdata.im         '); matching=0; end
if isequal(fdata.in          , data.output.in                   ) ; disp('ok: fdata.in         ') ; else disp('ERROR: fdata.in         '); matching=0; end
if isequal(fdata.lmns        , data.output.lmns                 ) ; disp('ok: fdata.lmns       ') ; else disp('ERROR: fdata.lmns       '); matching=0; end
if isequal(fdata.mn          , data.output.mn                   ) ; disp('ok: fdata.mn         ') ; else disp('ERROR: fdata.mn         '); matching=0; end
if isequal(fdata.volume      , data.output.volume               ) ; disp('ok: fdata.volume     ') ; else disp('ERROR: fdata.volume     '); matching=0; end

if isequal(fdata.Ate{1}      , data.vector_potential.Ate        ) ; disp('ok: fdata.Ate        ') ; else disp('ERROR: fdata.Ate        '); matching=0; end
if isequal(fdata.Aze{1}      , data.vector_potential.Aze        ) ; disp('ok: fdata.Aze        ') ; else disp('ERROR: fdata.Aze        '); matching=0; end
if isequal(fdata.Ato{1}      , data.vector_potential.Ato        ) ; disp('ok: fdata.Ato        ') ; else disp('ERROR: fdata.Ato        '); matching=0; end
if isequal(fdata.Azo{1}      , data.vector_potential.Azo        ) ; disp('ok: fdata.Azo        ') ; else disp('ERROR: fdata.Azo        '); matching=0; end

% check equivalence of gdata: from ".testcase.sp.grid" file
if isequal(gdata.Bnc         , data.output.Bnc                  ) ; disp('ok: gdata.Bnc        ') ; else disp('ERROR: gdata.Bnc        '); matching=0; end
if isequal(gdata.Bns         , data.output.Bns                  ) ; disp('ok: gdata.Bns        ') ; else disp('ERROR: gdata.Bns        '); matching=0; end
if isequal(gdata.Btemn       , data.output.Btemn                ) ; disp('ok: gdata.Btemn      ') ; else disp('ERROR: gdata.Btemn      '); matching=0; end
if isequal(gdata.Btomn       , data.output.Btomn                ) ; disp('ok: gdata.Btomn      ') ; else disp('ERROR: gdata.Btomn      '); matching=0; end
if isequal(gdata.Bzemn       , data.output.Bzemn                ) ; disp('ok: gdata.Bzemn      ') ; else disp('ERROR: gdata.Bzemn      '); matching=0; end
if isequal(gdata.Bzomn       , data.output.Bzomn                ) ; disp('ok: gdata.Bzomn      ') ; else disp('ERROR: gdata.Bzomn      '); matching=0; end
if isequal(gdata.ForceErr    , data.output.ForceErr             ) ; disp('ok: gdata.ForceErr   ') ; else disp('ERROR: gdata.ForceErr   '); matching=0; end
if isequal(gdata.Igeometry   , data.input.physics.Igeometry     ) ; disp('ok: gdata.Igeometry  ') ; else disp('ERROR: gdata.Igeometry  '); matching=0; end
if isequal(gdata.Istellsym   , data.input.physics.Istellsym     ) ; disp('ok: gdata.Istellsym  ') ; else disp('ERROR: gdata.Istellsym  '); matching=0; end
if isequal(gdata.Ladiabatic  , data.input.physics.Ladiabatic    ) ; disp('ok: gdata.Ladiabatic ') ; else disp('ERROR: gdata.Ladiabatic '); matching=0; end
if isequal(gdata.Lconstraint , data.input.physics.Lconstraint   ) ; disp('ok: gdata.Lconstraint') ; else disp('ERROR: gdata.Lconstraint'); matching=0; end
if isequal(gdata.Lfreebound  , data.input.physics.Lfreebound    ) ; disp('ok: gdata.Lfreebound ') ; else disp('ERROR: gdata.Lfreebound '); matching=0; end
if isequal(gdata.Lperturbed  , data.input.diagnostics.Lperturbed) ; disp('ok: gdata.Lperturbed ') ; else disp('ERROR: gdata.Lperturbed '); matching=0; end
if isequal(gdata.Lrad        , data.input.physics.Lrad          ) ; disp('ok: gdata.Lrad       ') ; else disp('ERROR: gdata.Lrad       '); matching=0; end
if isequal(gdata.Mpol        , data.input.physics.Mpol          ) ; disp('ok: gdata.Mpol       ') ; else disp('ERROR: gdata.Mpol       '); matching=0; end
if isequal(gdata.Mrad        , data.output.Mrad                 ) ; disp('ok: gdata.Mrad       ') ; else disp('ERROR: gdata.Mrad       '); matching=0; end
if isequal(gdata.Mvol        , data.output.Mvol                 ) ; disp('ok: gdata.Mvol       ') ; else disp('ERROR: gdata.Mvol       '); matching=0; end
if isequal(gdata.Nfp         , data.input.physics.Nfp           ) ; disp('ok: gdata.Nfp        ') ; else disp('ERROR: gdata.Nfp        '); matching=0; end
if isequal(gdata.Ntor        , data.input.physics.Ntor          ) ; disp('ok: gdata.Ntor       ') ; else disp('ERROR: gdata.Ntor       '); matching=0; end
if isequal(gdata.Nvol        , data.input.physics.Nvol          ) ; disp('ok: gdata.Nvol       ') ; else disp('ERROR: gdata.Nvol       '); matching=0; end
if isequal(gdata.Rbc         , data.output.Rbc                  ) ; disp('ok: gdata.Rbc        ') ; else disp('ERROR: gdata.Rbc        '); matching=0; end
if isequal(gdata.Rbs         , data.output.Rbs                  ) ; disp('ok: gdata.Rbs        ') ; else disp('ERROR: gdata.Rbs        '); matching=0; end
if isequal(gdata.TT          , data.output.TT                   ) ; disp('ok: gdata.TT         ') ; else disp('ERROR: gdata.TT         '); matching=0; end
if isequal(gdata.Vnc         , data.output.Vnc                  ) ; disp('ok: gdata.Vnc        ') ; else disp('ERROR: gdata.Vnc        '); matching=0; end
if isequal(gdata.Vns         , data.output.Vns                  ) ; disp('ok: gdata.Vns        ') ; else disp('ERROR: gdata.Vns        '); matching=0; end
if isequal(gdata.Zbc         , data.output.Zbc                  ) ; disp('ok: gdata.Zbc        ') ; else disp('ERROR: gdata.Zbc        '); matching=0; end
if isequal(gdata.Zbs         , data.output.Zbs                  ) ; disp('ok: gdata.Zbs        ') ; else disp('ERROR: gdata.Zbs        '); matching=0; end
if isequal(gdata.adiabatic   , data.output.adiabatic            ) ; disp('ok: gdata.adiabatic  ') ; else disp('ERROR: gdata.adiabatic  '); matching=0; end
if isequal(gdata.curpol      , data.input.physics.curpol        ) ; disp('ok: gdata.curpol     ') ; else disp('ERROR: gdata.curpol     '); matching=0; end
if isequal(gdata.curtor      , data.input.physics.curtor        ) ; disp('ok: gdata.curtor     ') ; else disp('ERROR: gdata.curtor     '); matching=0; end
if isequal(gdata.dpp         , data.input.diagnostics.dpp       ) ; disp('ok: gdata.dpp        ') ; else disp('ERROR: gdata.dpp        '); matching=0; end
if isequal(gdata.dqq         , data.input.diagnostics.dqq       ) ; disp('ok: gdata.dqq        ') ; else disp('ERROR: gdata.dqq        '); matching=0; end
if isequal(gdata.forcetol    , data.input.global.forcetol       ) ; disp('ok: gdata.forcetol   ') ; else disp('ERROR: gdata.forcetol   '); matching=0; end
if isequal(gdata.gamma       , data.input.physics.gamma         ) ; disp('ok: gdata.gamma      ') ; else disp('ERROR: gdata.gamma      '); matching=0; end
if isequal(gdata.helicity    , data.output.helicity             ) ; disp('ok: gdata.helicity   ') ; else disp('ERROR: gdata.helicity   '); matching=0; end
if isequal(gdata.im          , data.output.im                   ) ; disp('ok: gdata.im         ') ; else disp('ERROR: gdata.im         '); matching=0; end
if isequal(gdata.in          , data.output.in                   ) ; disp('ok: gdata.in         ') ; else disp('ERROR: gdata.in         '); matching=0; end
if isequal(gdata.iota        , data.input.physics.iota          ) ; disp('ok: gdata.iota       ') ; else disp('ERROR: gdata.iota       '); matching=0; end
if isequal(gdata.lmns        , data.output.lmns                 ) ; disp('ok: gdata.lmns       ') ; else disp('ERROR: gdata.lmns       '); matching=0; end
if isequal(gdata.lp          , data.input.physics.lp            ) ; disp('ok: gdata.lp         ') ; else disp('ERROR: gdata.lp         '); matching=0; end
if isequal(gdata.lq          , data.input.physics.lq            ) ; disp('ok: gdata.lq         ') ; else disp('ERROR: gdata.lq         '); matching=0; end
if isequal(gdata.mn          , data.output.mn                   ) ; disp('ok: gdata.mn         ') ; else disp('ERROR: gdata.mn         '); matching=0; end
if isequal(gdata.mu          , data.input.physics.mu            ) ; disp('ok: gdata.mu         ') ; else disp('ERROR: gdata.mu         '); matching=0; end
if isequal(gdata.mupfits     , data.input.physics.mupfits       ) ; disp('ok: gdata.mupfits    ') ; else disp('ERROR: gdata.mupfits    '); matching=0; end
if isequal(gdata.mupftol     , data.input.physics.mupftol       ) ; disp('ok: gdata.mupftol    ') ; else disp('ERROR: gdata.mupftol    '); matching=0; end
if isequal(gdata.oita        , data.input.physics.oita          ) ; disp('ok: gdata.oita       ') ; else disp('ERROR: gdata.oita       '); matching=0; end
if isequal(gdata.pflux       , data.input.physics.pflux         ) ; disp('ok: gdata.pflux      ') ; else disp('ERROR: gdata.pflux      '); matching=0; end
if isequal(gdata.phiedge     , data.input.physics.phiedge       ) ; disp('ok: gdata.phiedge    ') ; else disp('ERROR: gdata.phiedge    '); matching=0; end
if isequal(gdata.pl          , data.input.physics.pl            ) ; disp('ok: gdata.pl         ') ; else disp('ERROR: gdata.pl         '); matching=0; end
if isequal(gdata.pr          , data.input.physics.pr            ) ; disp('ok: gdata.pr         ') ; else disp('ERROR: gdata.pr         '); matching=0; end
if isequal(gdata.pressure    , data.input.physics.pressure      ) ; disp('ok: gdata.pressure   ') ; else disp('ERROR: gdata.pressure   '); matching=0; end
if isequal(gdata.pscale      , data.input.physics.pscale        ) ; disp('ok: gdata.pscale     ') ; else disp('ERROR: gdata.pscale     '); matching=0; end
if isequal(gdata.ql          , data.input.physics.ql            ) ; disp('ok: gdata.ql         ') ; else disp('ERROR: gdata.ql         '); matching=0; end
if isequal(gdata.qr          , data.input.physics.qr            ) ; disp('ok: gdata.qr         ') ; else disp('ERROR: gdata.qr         '); matching=0; end
if isequal(gdata.rp          , data.input.physics.rp            ) ; disp('ok: gdata.rp         ') ; else disp('ERROR: gdata.rp         '); matching=0; end
if isequal(gdata.rq          , data.input.physics.rq            ) ; disp('ok: gdata.rq         ') ; else disp('ERROR: gdata.rq         '); matching=0; end
if isequal(gdata.tflux       , data.input.physics.tflux         ) ; disp('ok: gdata.tflux      ') ; else disp('ERROR: gdata.tflux      '); matching=0; end
if isequal(gdata.volume      , data.output.volume               ) ; disp('ok: gdata.volume     ') ; else disp('ERROR: gdata.volume     '); matching=0; end

if isequal(gdata.Nt          , data.grid.Nt                     ) ; disp('ok: gdata.Nt         ') ; else disp('ERROR: gdata.Nt         '); matching=0; end
if isequal(gdata.Nz          , data.grid.Nz                     ) ; disp('ok: gdata.Nz         ') ; else disp('ERROR: gdata.Nz         '); matching=0; end
if isequal(gdata.Ntz         , data.grid.Ntz                    ) ; disp('ok: gdata.Ntz        ') ; else disp('ERROR: gdata.Ntz        '); matching=0; end
if isequal(squeeze(gdata.Rij), data.grid.Rij'                   ) ; disp('ok: gdata.Rij        ') ; else disp('ERROR: gdata.Rij        '); matching=0; end
if isequal(squeeze(gdata.Zij), data.grid.Zij'                   ) ; disp('ok: gdata.Zij        ') ; else disp('ERROR: gdata.Zij        '); matching=0; end
if isequal(squeeze(gdata.sg) , data.grid.sg'                    ) ; disp('ok: gdata.sg         ') ; else disp('ERROR: gdata.sg         '); matching=0; end
if isequal(squeeze(gdata.BR) , data.grid.BR'                    ) ; disp('ok: gdata.BR         ') ; else disp('ERROR: gdata.BR         '); matching=0; end
if isequal(squeeze(gdata.Bp) , data.grid.Bp'                    ) ; disp('ok: gdata.Bp         ') ; else disp('ERROR: gdata.Bp         '); matching=0; end
if isequal(squeeze(gdata.BZ) , data.grid.BZ'                    ) ; disp('ok: gdata.BZ         ') ; else disp('ERROR: gdata.BZ         '); matching=0; end

% data on rotational transform from .ext.sp.t.0001.dat and so on
if isequal(idata.Bnc         , data.output.Bnc                  ) ; disp('ok: idata.Bnc        ') ; else disp('ERROR: idata.Bnc        '); matching=0; end
if isequal(idata.Bns         , data.output.Bns                  ) ; disp('ok: idata.Bns        ') ; else disp('ERROR: idata.Bns        '); matching=0; end
if isequal(idata.Btemn       , data.output.Btemn                ) ; disp('ok: idata.Btemn      ') ; else disp('ERROR: idata.Btemn      '); matching=0; end
if isequal(idata.Btomn       , data.output.Btomn                ) ; disp('ok: idata.Btomn      ') ; else disp('ERROR: idata.Btomn      '); matching=0; end
if isequal(idata.Bzemn       , data.output.Bzemn                ) ; disp('ok: idata.Bzemn      ') ; else disp('ERROR: idata.Bzemn      '); matching=0; end
if isequal(idata.Bzomn       , data.output.Bzomn                ) ; disp('ok: idata.Bzomn      ') ; else disp('ERROR: idata.Bzomn      '); matching=0; end
if isequal(idata.ForceErr    , data.output.ForceErr             ) ; disp('ok: idata.ForceErr   ') ; else disp('ERROR: idata.ForceErr   '); matching=0; end
if isequal(idata.Igeometry   , data.input.physics.Igeometry     ) ; disp('ok: idata.Igeometry  ') ; else disp('ERROR: idata.Igeometry  '); matching=0; end
if isequal(idata.Istellsym   , data.input.physics.Istellsym     ) ; disp('ok: idata.Istellsym  ') ; else disp('ERROR: idata.Istellsym  '); matching=0; end
if isequal(idata.Ladiabatic  , data.input.physics.Ladiabatic    ) ; disp('ok: idata.Ladiabatic ') ; else disp('ERROR: idata.Ladiabatic '); matching=0; end
if isequal(idata.Lconstraint , data.input.physics.Lconstraint   ) ; disp('ok: idata.Lconstraint') ; else disp('ERROR: idata.Lconstraint'); matching=0; end
if isequal(idata.Lfreebound  , data.input.physics.Lfreebound    ) ; disp('ok: idata.Lfreebound ') ; else disp('ERROR: idata.Lfreebound '); matching=0; end
if isequal(idata.Lperturbed  , data.input.diagnostics.Lperturbed) ; disp('ok: idata.Lperturbed ') ; else disp('ERROR: idata.Lperturbed '); matching=0; end
if isequal(idata.Lrad        , data.input.physics.Lrad          ) ; disp('ok: idata.Lrad       ') ; else disp('ERROR: idata.Lrad       '); matching=0; end
if isequal(idata.Mpol        , data.input.physics.Mpol          ) ; disp('ok: idata.Mpol       ') ; else disp('ERROR: idata.Mpol       '); matching=0; end
if isequal(idata.Mrad        , data.output.Mrad                 ) ; disp('ok: idata.Mrad       ') ; else disp('ERROR: idata.Mrad       '); matching=0; end
if isequal(idata.Mvol        , data.output.Mvol                 ) ; disp('ok: idata.Mvol       ') ; else disp('ERROR: idata.Mvol       '); matching=0; end
if isequal(idata.Nfp         , data.input.physics.Nfp           ) ; disp('ok: idata.Nfp        ') ; else disp('ERROR: idata.Nfp        '); matching=0; end
if isequal(idata.Ntor        , data.input.physics.Ntor          ) ; disp('ok: idata.Ntor       ') ; else disp('ERROR: idata.Ntor       '); matching=0; end
if isequal(idata.Nvol        , data.input.physics.Nvol          ) ; disp('ok: idata.Nvol       ') ; else disp('ERROR: idata.Nvol       '); matching=0; end
if isequal(idata.Rbc         , data.output.Rbc                  ) ; disp('ok: idata.Rbc        ') ; else disp('ERROR: idata.Rbc        '); matching=0; end
if isequal(idata.Rbs         , data.output.Rbs                  ) ; disp('ok: idata.Rbs        ') ; else disp('ERROR: idata.Rbs        '); matching=0; end
if isequal(idata.TT          , data.output.TT                   ) ; disp('ok: idata.TT         ') ; else disp('ERROR: idata.TT         '); matching=0; end
if isequal(idata.Vnc         , data.output.Vnc                  ) ; disp('ok: idata.Vnc        ') ; else disp('ERROR: idata.Vnc        '); matching=0; end
if isequal(idata.Vns         , data.output.Vns                  ) ; disp('ok: idata.Vns        ') ; else disp('ERROR: idata.Vns        '); matching=0; end
if isequal(idata.Zbc         , data.output.Zbc                  ) ; disp('ok: idata.Zbc        ') ; else disp('ERROR: idata.Zbc        '); matching=0; end
if isequal(idata.Zbs         , data.output.Zbs                  ) ; disp('ok: idata.Zbs        ') ; else disp('ERROR: idata.Zbs        '); matching=0; end
if isequal(idata.adiabatic   , data.output.adiabatic            ) ; disp('ok: idata.adiabatic  ') ; else disp('ERROR: idata.adiabatic  '); matching=0; end
if isequal(idata.curpol      , data.input.physics.curpol        ) ; disp('ok: idata.curpol     ') ; else disp('ERROR: idata.curpol     '); matching=0; end
if isequal(idata.curtor      , data.input.physics.curtor        ) ; disp('ok: idata.curtor     ') ; else disp('ERROR: idata.curtor     '); matching=0; end
if isequal(idata.dpp         , data.input.diagnostics.dpp       ) ; disp('ok: idata.dpp        ') ; else disp('ERROR: idata.dpp        '); matching=0; end
if isequal(idata.dqq         , data.input.diagnostics.dqq       ) ; disp('ok: idata.dqq        ') ; else disp('ERROR: idata.dqq        '); matching=0; end
if isequal(idata.forcetol    , data.input.global.forcetol       ) ; disp('ok: idata.forcetol   ') ; else disp('ERROR: idata.forcetol   '); matching=0; end
if isequal(idata.gamma       , data.input.physics.gamma         ) ; disp('ok: idata.gamma      ') ; else disp('ERROR: idata.gamma      '); matching=0; end
if isequal(idata.helicity    , data.output.helicity             ) ; disp('ok: idata.helicity   ') ; else disp('ERROR: idata.helicity   '); matching=0; end
if isequal(idata.im          , data.output.im                   ) ; disp('ok: idata.im         ') ; else disp('ERROR: idata.im         '); matching=0; end
if isequal(idata.in          , data.output.in                   ) ; disp('ok: idata.in         ') ; else disp('ERROR: idata.in         '); matching=0; end
if isequal(idata.iota        , data.input.physics.iota          ) ; disp('ok: idata.iota       ') ; else disp('ERROR: idata.iota       '); matching=0; end
if isequal(idata.lmns        , data.output.lmns                 ) ; disp('ok: idata.lmns       ') ; else disp('ERROR: idata.lmns       '); matching=0; end
if isequal(idata.lp          , data.input.physics.lp            ) ; disp('ok: idata.lp         ') ; else disp('ERROR: idata.lp         '); matching=0; end
if isequal(idata.lq          , data.input.physics.lq            ) ; disp('ok: idata.lq         ') ; else disp('ERROR: idata.lq         '); matching=0; end
if isequal(idata.mn          , data.output.mn                   ) ; disp('ok: idata.mn         ') ; else disp('ERROR: idata.mn         '); matching=0; end
if isequal(idata.mu          , data.input.physics.mu            ) ; disp('ok: idata.mu         ') ; else disp('ERROR: idata.mu         '); matching=0; end
if isequal(idata.mupfits     , data.input.physics.mupfits       ) ; disp('ok: idata.mupfits    ') ; else disp('ERROR: idata.mupfits    '); matching=0; end
if isequal(idata.mupftol     , data.input.physics.mupftol       ) ; disp('ok: idata.mupftol    ') ; else disp('ERROR: idata.mupftol    '); matching=0; end
if isequal(idata.oita        , data.input.physics.oita          ) ; disp('ok: idata.oita       ') ; else disp('ERROR: idata.oita       '); matching=0; end
if isequal(idata.pflux       , data.input.physics.pflux         ) ; disp('ok: idata.pflux      ') ; else disp('ERROR: idata.pflux      '); matching=0; end
if isequal(idata.phiedge     , data.input.physics.phiedge       ) ; disp('ok: idata.phiedge    ') ; else disp('ERROR: idata.phiedge    '); matching=0; end
if isequal(idata.pl          , data.input.physics.pl            ) ; disp('ok: idata.pl         ') ; else disp('ERROR: idata.pl         '); matching=0; end
if isequal(idata.pr          , data.input.physics.pr            ) ; disp('ok: idata.pr         ') ; else disp('ERROR: idata.pr         '); matching=0; end
if isequal(idata.pressure    , data.input.physics.pressure      ) ; disp('ok: idata.pressure   ') ; else disp('ERROR: idata.pressure   '); matching=0; end
if isequal(idata.pscale      , data.input.physics.pscale        ) ; disp('ok: idata.pscale     ') ; else disp('ERROR: idata.pscale     '); matching=0; end
if isequal(idata.ql          , data.input.physics.ql            ) ; disp('ok: idata.ql         ') ; else disp('ERROR: idata.ql         '); matching=0; end
if isequal(idata.qr          , data.input.physics.qr            ) ; disp('ok: idata.qr         ') ; else disp('ERROR: idata.qr         '); matching=0; end
if isequal(idata.rp          , data.input.physics.rp            ) ; disp('ok: idata.rp         ') ; else disp('ERROR: idata.rp         '); matching=0; end
if isequal(idata.rq          , data.input.physics.rq            ) ; disp('ok: idata.rq         ') ; else disp('ERROR: idata.rq         '); matching=0; end
if isequal(idata.tflux       , data.input.physics.tflux         ) ; disp('ok: idata.tflux      ') ; else disp('ERROR: idata.tflux      '); matching=0; end
if isequal(idata.volume      , data.output.volume               ) ; disp('ok: idata.volume     ') ; else disp('ERROR: idata.volume     '); matching=0; end

if isequal(idata.sarr        , data.transform.fiota(:,1)'       ) ; disp('ok: idata.sarr       ') ; else disp('ERROR: idata.sarr       '); matching=0; end





if (matching == 0)
  disp('Not maching :(')
else
  disp('Matching :)')
end

end

