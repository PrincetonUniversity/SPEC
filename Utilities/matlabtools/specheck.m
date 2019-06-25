<<<<<<< HEAD

=======
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

% quantities from old .sp.h5 file, sorted by equivalent quantity in the new output file
if isequal(fdata.forcetol   , data.input.global.forcetol       ) ; disp('ok: fdata.forcetol   ') ; else disp('ERROR: fdata.forcetol   '); matching=0; end

if isequal(fdata.Igeometry  , data.input.physics.Igeometry     ) ; disp('ok: fdata.Igeometry  ') ; else disp('ERROR: fdata.Igeometry  '); matching=0; end
if isequal(fdata.Istellsym  , data.input.physics.Istellsym     ) ; disp('ok: fdata.Istellsym  ') ; else disp('ERROR: fdata.Istellsym  '); matching=0; end
if isequal(fdata.Ladiabatic , data.input.physics.Ladiabatic    ) ; disp('ok: fdata.Ladiabatic ') ; else disp('ERROR: fdata.Ladiabatic '); matching=0; end
if isequal(fdata.Lconstraint, data.input.physics.Lconstraint   ) ; disp('ok: fdata.Lconstraint') ; else disp('ERROR: fdata.Lconstraint'); matching=0; end
if isequal(fdata.Lfreebound , data.input.physics.Lfreebound    ) ; disp('ok: fdata.Lfreebound ') ; else disp('ERROR: fdata.Lfreebound '); matching=0; end
if isequal(fdata.Lrad       , data.input.physics.Lrad          ) ; disp('ok: fdata.Lrad       ') ; else disp('ERROR: fdata.Lrad       '); matching=0; end
if isequal(fdata.Mpol       , data.input.physics.Mpol          ) ; disp('ok: fdata.Mpol       ') ; else disp('ERROR: fdata.Mpol       '); matching=0; end
if isequal(fdata.Nfp        , data.input.physics.Nfp           ) ; disp('ok: fdata.Nfp        ') ; else disp('ERROR: fdata.Nfp        '); matching=0; end
if isequal(fdata.Ntor       , data.input.physics.Ntor          ) ; disp('ok: fdata.Ntor       ') ; else disp('ERROR: fdata.Ntor       '); matching=0; end
if isequal(fdata.Nvol       , data.input.physics.Nvol          ) ; disp('ok: fdata.Nvol       ') ; else disp('ERROR: fdata.Nvol       '); matching=0; end
if isequal(fdata.curpol     , data.input.physics.curpol        ) ; disp('ok: fdata.curpol     ') ; else disp('ERROR: fdata.curpol     '); matching=0; end
if isequal(fdata.curtor     , data.input.physics.curtor        ) ; disp('ok: fdata.curtor     ') ; else disp('ERROR: fdata.curtor     '); matching=0; end
if isequal(fdata.gamma      , data.input.physics.gamma         ) ; disp('ok: fdata.gamma      ') ; else disp('ERROR: fdata.gamma      '); matching=0; end
if isequal(fdata.iota       , data.input.physics.iota          ) ; disp('ok: fdata.iota       ') ; else disp('ERROR: fdata.iota       '); matching=0; end
if isequal(fdata.lp         , data.input.physics.lp            ) ; disp('ok: fdata.lp         ') ; else disp('ERROR: fdata.lp         '); matching=0; end
if isequal(fdata.lq         , data.input.physics.lq            ) ; disp('ok: fdata.lq         ') ; else disp('ERROR: fdata.lq         '); matching=0; end
if isequal(fdata.mu         , data.input.physics.mu            ) ; disp('ok: fdata.mu         ') ; else disp('ERROR: fdata.mu         '); matching=0; end
if isequal(fdata.mupfits    , data.input.physics.mupfits       ) ; disp('ok: fdata.mupfits    ') ; else disp('ERROR: fdata.mupfits    '); matching=0; end
if isequal(fdata.mupftol    , data.input.physics.mupftol       ) ; disp('ok: fdata.mupftol    ') ; else disp('ERROR: fdata.mupftol    '); matching=0; end
if isequal(fdata.oita       , data.input.physics.oita          ) ; disp('ok: fdata.oita       ') ; else disp('ERROR: fdata.oita       '); matching=0; end
if isequal(fdata.pflux      , data.input.physics.pflux         ) ; disp('ok: fdata.pflux      ') ; else disp('ERROR: fdata.pflux      '); matching=0; end
if isequal(fdata.phiedge    , data.input.physics.phiedge       ) ; disp('ok: fdata.phiedge    ') ; else disp('ERROR: fdata.phiedge    '); matching=0; end
if isequal(fdata.pl         , data.input.physics.pl            ) ; disp('ok: fdata.pl         ') ; else disp('ERROR: fdata.pl         '); matching=0; end
if isequal(fdata.pr         , data.input.physics.pr            ) ; disp('ok: fdata.pr         ') ; else disp('ERROR: fdata.pr         '); matching=0; end
if isequal(fdata.pressure   , data.input.physics.pressure      ) ; disp('ok: fdata.pressure   ') ; else disp('ERROR: fdata.pressure   '); matching=0; end
if isequal(fdata.pscale     , data.input.physics.pscale        ) ; disp('ok: fdata.pscale     ') ; else disp('ERROR: fdata.pscale     '); matching=0; end
if isequal(fdata.ql         , data.input.physics.ql            ) ; disp('ok: fdata.ql         ') ; else disp('ERROR: fdata.ql         '); matching=0; end
if isequal(fdata.qr         , data.input.physics.qr            ) ; disp('ok: fdata.qr         ') ; else disp('ERROR: fdata.qr         '); matching=0; end
if isequal(fdata.rp         , data.input.physics.rp            ) ; disp('ok: fdata.rp         ') ; else disp('ERROR: fdata.rp         '); matching=0; end
if isequal(fdata.rq         , data.input.physics.rq            ) ; disp('ok: fdata.rq         ') ; else disp('ERROR: fdata.rq         '); matching=0; end
if isequal(fdata.tflux      , data.input.physics.tflux         ) ; disp('ok: fdata.tflux      ') ; else disp('ERROR: fdata.tflux      '); matching=0; end

if isequal(fdata.Lperturbed , data.input.diagnostics.Lperturbed) ; disp('ok: fdata.Lperturbed ') ; else disp('ERROR: fdata.Lperturbed '); matching=0; end
if isequal(fdata.dpp        , data.input.diagnostics.dpp       ) ; disp('ok: fdata.dpp        ') ; else disp('ERROR: fdata.dpp        '); matching=0; end
if isequal(fdata.dqq        , data.input.diagnostics.dqq       ) ; disp('ok: fdata.dqq        ') ; else disp('ERROR: fdata.dqq        '); matching=0; end

if isequal(fdata.Bnc        , data.output.Bnc                  ) ; disp('ok: fdata.Bnc        ') ; else disp('ERROR: fdata.Bnc        '); matching=0; end
if isequal(fdata.Bns        , data.output.Bns                  ) ; disp('ok: fdata.Bns        ') ; else disp('ERROR: fdata.Bns        '); matching=0; end
if isequal(fdata.Btemn      , data.output.Btemn                ) ; disp('ok: fdata.Btemn      ') ; else disp('ERROR: fdata.Btemn      '); matching=0; end
if isequal(fdata.Btomn      , data.output.Btomn                ) ; disp('ok: fdata.Btomn      ') ; else disp('ERROR: fdata.Btomn      '); matching=0; end
if isequal(fdata.Bzemn      , data.output.Bzemn                ) ; disp('ok: fdata.Bzemn      ') ; else disp('ERROR: fdata.Bzemn      '); matching=0; end
if isequal(fdata.Bzomn      , data.output.Bzomn                ) ; disp('ok: fdata.Bzomn      ') ; else disp('ERROR: fdata.Bzomn      '); matching=0; end
if isequal(fdata.ForceErr   , data.output.ForceErr             ) ; disp('ok: fdata.ForceErr   ') ; else disp('ERROR: fdata.ForceErr   '); matching=0; end
if isequal(fdata.Mrad       , data.output.Mrad                 ) ; disp('ok: fdata.Mrad       ') ; else disp('ERROR: fdata.Mrad       '); matching=0; end
if isequal(fdata.Mvol       , data.output.Mvol                 ) ; disp('ok: fdata.Mvol       ') ; else disp('ERROR: fdata.Mvol       '); matching=0; end
if isequal(fdata.Rbc        , data.output.Rbc                  ) ; disp('ok: fdata.Rbc        ') ; else disp('ERROR: fdata.Rbc        '); matching=0; end
if isequal(fdata.Rbs        , data.output.Rbs                  ) ; disp('ok: fdata.Rbs        ') ; else disp('ERROR: fdata.Rbs        '); matching=0; end
if isequal(fdata.TT         , data.output.TT                   ) ; disp('ok: fdata.TT         ') ; else disp('ERROR: fdata.TT         '); matching=0; end
if isequal(fdata.Vnc        , data.output.Vnc                  ) ; disp('ok: fdata.Vnc        ') ; else disp('ERROR: fdata.Vnc        '); matching=0; end
if isequal(fdata.Vns        , data.output.Vns                  ) ; disp('ok: fdata.Vns        ') ; else disp('ERROR: fdata.Vns        '); matching=0; end
if isequal(fdata.Zbc        , data.output.Zbc                  ) ; disp('ok: fdata.Zbc        ') ; else disp('ERROR: fdata.Zbc        '); matching=0; end
if isequal(fdata.Zbs        , data.output.Zbs                  ) ; disp('ok: fdata.Zbs        ') ; else disp('ERROR: fdata.Zbs        '); matching=0; end
if isequal(fdata.adiabatic  , data.output.adiabatic            ) ; disp('ok: fdata.adiabatic  ') ; else disp('ERROR: fdata.adiabatic  '); matching=0; end
if isequal(fdata.helicity   , data.output.helicity             ) ; disp('ok: fdata.helicity   ') ; else disp('ERROR: fdata.helicity   '); matching=0; end
if isequal(fdata.im         , data.output.im                   ) ; disp('ok: fdata.im         ') ; else disp('ERROR: fdata.im         '); matching=0; end
if isequal(fdata.in         , data.output.in                   ) ; disp('ok: fdata.in         ') ; else disp('ERROR: fdata.in         '); matching=0; end
if isequal(fdata.lmns       , data.output.lmns                 ) ; disp('ok: fdata.lmns       ') ; else disp('ERROR: fdata.lmns       '); matching=0; end
if isequal(fdata.mn         , data.output.mn                   ) ; disp('ok: fdata.mn         ') ; else disp('ERROR: fdata.mn         '); matching=0; end
if isequal(fdata.volume     , data.output.volume               ) ; disp('ok: fdata.volume     ') ; else disp('ERROR: fdata.volume     '); matching=0; end

if isequal(fdata.Ate{1}     , data.vector_potential.Ate        ) ; disp('ok: fdata.Ate        ') ; else disp('ERROR: fdata.Ate        '); matching=0; end
if isequal(fdata.Aze{1}     , data.vector_potential.Aze        ) ; disp('ok: fdata.Aze        ') ; else disp('ERROR: fdata.Aze        '); matching=0; end
if isequal(fdata.Ato{1}     , data.vector_potential.Ato        ) ; disp('ok: fdata.Ato        ') ; else disp('ERROR: fdata.Ato        '); matching=0; end
if isequal(fdata.Azo{1}     , data.vector_potential.Azo        ) ; disp('ok: fdata.Azo        ') ; else disp('ERROR: fdata.Azo        '); matching=0; end

if (matching == 0)
  disp('Not maching :(')
else
  disp('Matching :)')
end

end

