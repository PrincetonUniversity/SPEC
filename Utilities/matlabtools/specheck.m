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

tol=1e-12;

% check equivalence of fdata: from "testcase.sp.h5" file
% quantities from old .sp.h5 file, sorted by equivalent quantity in the new output file
if isequaltol(fdata.forcetol    , data.input.global.forcetol       , tol) ; disp('ok: fdata.forcetol   ') ; else ; disp('ERROR: fdata.forcetol   '); matching=0; end

if isequaltol(fdata.Igeometry   , data.input.physics.Igeometry     , tol) ; disp('ok: fdata.Igeometry  ') ; else ; disp('ERROR: fdata.Igeometry  '); matching=0; end
if isequaltol(fdata.Istellsym   , data.input.physics.Istellsym     , tol) ; disp('ok: fdata.Istellsym  ') ; else ; disp('ERROR: fdata.Istellsym  '); matching=0; end
if isequaltol(fdata.Ladiabatic  , data.input.physics.Ladiabatic    , tol) ; disp('ok: fdata.Ladiabatic ') ; else ; disp('ERROR: fdata.Ladiabatic '); matching=0; end
if isequaltol(fdata.Lconstraint , data.input.physics.Lconstraint   , tol) ; disp('ok: fdata.Lconstraint') ; else ; disp('ERROR: fdata.Lconstraint'); matching=0; end
if isequaltol(fdata.Lfreebound  , data.input.physics.Lfreebound    , tol) ; disp('ok: fdata.Lfreebound ') ; else ; disp('ERROR: fdata.Lfreebound '); matching=0; end
if isequaltol(fdata.Lrad        , data.input.physics.Lrad          , tol) ; disp('ok: fdata.Lrad       ') ; else ; disp('ERROR: fdata.Lrad       '); matching=0; end
if isequaltol(fdata.Mpol        , data.input.physics.Mpol          , tol) ; disp('ok: fdata.Mpol       ') ; else ; disp('ERROR: fdata.Mpol       '); matching=0; end
if isequaltol(fdata.Nfp         , data.input.physics.Nfp           , tol) ; disp('ok: fdata.Nfp        ') ; else ; disp('ERROR: fdata.Nfp        '); matching=0; end
if isequaltol(fdata.Ntor        , data.input.physics.Ntor          , tol) ; disp('ok: fdata.Ntor       ') ; else ; disp('ERROR: fdata.Ntor       '); matching=0; end
if isequaltol(fdata.Nvol        , data.input.physics.Nvol          , tol) ; disp('ok: fdata.Nvol       ') ; else ; disp('ERROR: fdata.Nvol       '); matching=0; end
if isequaltol(fdata.curpol      , data.input.physics.curpol        , tol) ; disp('ok: fdata.curpol     ') ; else ; disp('ERROR: fdata.curpol     '); matching=0; end
if isequaltol(fdata.curtor      , data.input.physics.curtor        , tol) ; disp('ok: fdata.curtor     ') ; else ; disp('ERROR: fdata.curtor     '); matching=0; end
if isequaltol(fdata.gamma       , data.input.physics.gamma         , tol) ; disp('ok: fdata.gamma      ') ; else ; disp('ERROR: fdata.gamma      '); matching=0; end
if isequaltol(fdata.iota        , data.input.physics.iota          , tol) ; disp('ok: fdata.iota       ') ; else ; disp('ERROR: fdata.iota       '); matching=0; end
if isequaltol(fdata.lp          , data.input.physics.lp            , tol) ; disp('ok: fdata.lp         ') ; else ; disp('ERROR: fdata.lp         '); matching=0; end
if isequaltol(fdata.lq          , data.input.physics.lq            , tol) ; disp('ok: fdata.lq         ') ; else ; disp('ERROR: fdata.lq         '); matching=0; end
if isequaltol(fdata.mupfits     , data.input.physics.mupfits       , tol) ; disp('ok: fdata.mupfits    ') ; else ; disp('ERROR: fdata.mupfits    '); matching=0; end
if isequaltol(fdata.mupftol     , data.input.physics.mupftol       , tol) ; disp('ok: fdata.mupftol    ') ; else ; disp('ERROR: fdata.mupftol    '); matching=0; end
if isequaltol(fdata.oita        , data.input.physics.oita          , tol) ; disp('ok: fdata.oita       ') ; else ; disp('ERROR: fdata.oita       '); matching=0; end
if isequaltol(fdata.phiedge     , data.input.physics.phiedge       , tol) ; disp('ok: fdata.phiedge    ') ; else ; disp('ERROR: fdata.phiedge    '); matching=0; end
if isequaltol(fdata.pl          , data.input.physics.pl            , tol) ; disp('ok: fdata.pl         ') ; else ; disp('ERROR: fdata.pl         '); matching=0; end
if isequaltol(fdata.pr          , data.input.physics.pr            , tol) ; disp('ok: fdata.pr         ') ; else ; disp('ERROR: fdata.pr         '); matching=0; end
if isequaltol(fdata.pressure    , data.input.physics.pressure      , tol) ; disp('ok: fdata.pressure   ') ; else ; disp('ERROR: fdata.pressure   '); matching=0; end
if isequaltol(fdata.pscale      , data.input.physics.pscale        , tol) ; disp('ok: fdata.pscale     ') ; else ; disp('ERROR: fdata.pscale     '); matching=0; end
if isequaltol(fdata.ql          , data.input.physics.ql            , tol) ; disp('ok: fdata.ql         ') ; else ; disp('ERROR: fdata.ql         '); matching=0; end
if isequaltol(fdata.qr          , data.input.physics.qr            , tol) ; disp('ok: fdata.qr         ') ; else ; disp('ERROR: fdata.qr         '); matching=0; end
if isequaltol(fdata.rp          , data.input.physics.rp            , tol) ; disp('ok: fdata.rp         ') ; else ; disp('ERROR: fdata.rp         '); matching=0; end
if isequaltol(fdata.rq          , data.input.physics.rq            , tol) ; disp('ok: fdata.rq         ') ; else ; disp('ERROR: fdata.rq         '); matching=0; end

if isequaltol(fdata.Lperturbed  , data.input.diagnostics.Lperturbed, tol) ; disp('ok: fdata.Lperturbed ') ; else ; disp('ERROR: fdata.Lperturbed '); matching=0; end
if isequaltol(fdata.dpp         , data.input.diagnostics.dpp       , tol) ; disp('ok: fdata.dpp        ') ; else ; disp('ERROR: fdata.dpp        '); matching=0; end
if isequaltol(fdata.dqq         , data.input.diagnostics.dqq       , tol) ; disp('ok: fdata.dqq        ') ; else ; disp('ERROR: fdata.dqq        '); matching=0; end

if isequaltol(fdata.Bnc         , data.output.Bnc                  , tol) ; disp('ok: fdata.Bnc        ') ; else ; disp('ERROR: fdata.Bnc        '); matching=0; end
if isequaltol(fdata.Bns         , data.output.Bns                  , tol) ; disp('ok: fdata.Bns        ') ; else ; disp('ERROR: fdata.Bns        '); matching=0; end
if isequaltol(fdata.Btemn       , data.output.Btemn                , tol) ; disp('ok: fdata.Btemn      ') ; else ; disp('ERROR: fdata.Btemn      '); matching=0; end
if isequaltol(fdata.Btomn       , data.output.Btomn                , tol) ; disp('ok: fdata.Btomn      ') ; else ; disp('ERROR: fdata.Btomn      '); matching=0; end
if isequaltol(fdata.Bzemn       , data.output.Bzemn                , tol) ; disp('ok: fdata.Bzemn      ') ; else ; disp('ERROR: fdata.Bzemn      '); matching=0; end
if isequaltol(fdata.Bzomn       , data.output.Bzomn                , tol) ; disp('ok: fdata.Bzomn      ') ; else ; disp('ERROR: fdata.Bzomn      '); matching=0; end
if isequaltol(fdata.ForceErr    , data.output.ForceErr             , tol) ; disp('ok: fdata.ForceErr   ') ; else ; disp('ERROR: fdata.ForceErr   '); matching=0; end
if isequaltol(fdata.Mrad        , data.output.Mrad                 , tol) ; disp('ok: fdata.Mrad       ') ; else ; disp('ERROR: fdata.Mrad       '); matching=0; end
if isequaltol(fdata.Mvol        , data.output.Mvol                 , tol) ; disp('ok: fdata.Mvol       ') ; else ; disp('ERROR: fdata.Mvol       '); matching=0; end
if isequaltol(fdata.Rbc         , data.output.Rbc                  , tol) ; disp('ok: fdata.Rbc        ') ; else ; disp('ERROR: fdata.Rbc        '); matching=0; end
if isequaltol(fdata.Rbs         , data.output.Rbs                  , tol) ; disp('ok: fdata.Rbs        ') ; else ; disp('ERROR: fdata.Rbs        '); matching=0; end
if isequaltol(fdata.TT          , data.output.TT                   , tol) ; disp('ok: fdata.TT         ') ; else ; disp('ERROR: fdata.TT         '); matching=0; end
if isequaltol(fdata.Vnc         , data.output.Vnc                  , tol) ; disp('ok: fdata.Vnc        ') ; else ; disp('ERROR: fdata.Vnc        '); matching=0; end
if isequaltol(fdata.Vns         , data.output.Vns                  , tol) ; disp('ok: fdata.Vns        ') ; else ; disp('ERROR: fdata.Vns        '); matching=0; end
if isequaltol(fdata.Zbc         , data.output.Zbc                  , tol) ; disp('ok: fdata.Zbc        ') ; else ; disp('ERROR: fdata.Zbc        '); matching=0; end
if isequaltol(fdata.Zbs         , data.output.Zbs                  , tol) ; disp('ok: fdata.Zbs        ') ; else ; disp('ERROR: fdata.Zbs        '); matching=0; end
if isequaltol(fdata.adiabatic   , data.output.adiabatic            , tol) ; disp('ok: fdata.adiabatic  ') ; else ; disp('ERROR: fdata.adiabatic  '); matching=0; end
if isequaltol(fdata.helicity    , data.output.helicity             , tol) ; disp('ok: fdata.helicity   ') ; else ; disp('ERROR: fdata.helicity   '); matching=0; end
if isequaltol(fdata.im          , data.output.im                   , tol) ; disp('ok: fdata.im         ') ; else ; disp('ERROR: fdata.im         '); matching=0; end
if isequaltol(fdata.in          , data.output.in                   , tol) ; disp('ok: fdata.in         ') ; else ; disp('ERROR: fdata.in         '); matching=0; end
if isequaltol(fdata.lmns        , data.output.lmns                 , tol) ; disp('ok: fdata.lmns       ') ; else ; disp('ERROR: fdata.lmns       '); matching=0; end
if isequaltol(fdata.mn          , data.output.mn                   , tol) ; disp('ok: fdata.mn         ') ; else ; disp('ERROR: fdata.mn         '); matching=0; end
if isequaltol(fdata.volume      , data.output.volume               , tol) ; disp('ok: fdata.volume     ') ; else ; disp('ERROR: fdata.volume     '); matching=0; end
if isequaltol(fdata.mu          , data.output.mu                   , tol) ; disp('ok: fdata.mu         ') ; else ; disp('ERROR: fdata.mu         '); matching=0; end
if isequaltol(fdata.pflux       , data.output.pflux                , tol) ; disp('ok: fdata.pflux      ') ; else ; disp('ERROR: fdata.pflux      '); matching=0; end
if isequaltol(fdata.tflux       , data.output.tflux                , tol) ; disp('ok: fdata.tflux      ') ; else ; disp('ERROR: fdata.tflux      '); matching=0; end

if isequaltol(fdata.Ate         , data.vector_potential.Ate        , tol) ; disp('ok: fdata.Ate        ') ; else ; disp('ERROR: fdata.Ate        '); matching=0; end
if isequaltol(fdata.Aze         , data.vector_potential.Aze        , tol) ; disp('ok: fdata.Aze        ') ; else ; disp('ERROR: fdata.Aze        '); matching=0; end
if isequaltol(fdata.Ato         , data.vector_potential.Ato        , tol) ; disp('ok: fdata.Ato        ') ; else ; disp('ERROR: fdata.Ato        '); matching=0; end
if isequaltol(fdata.Azo         , data.vector_potential.Azo        , tol) ; disp('ok: fdata.Azo        ') ; else ; disp('ERROR: fdata.Azo        '); matching=0; end

% check equivalence of gdata: from ".testcase.sp.grid" file
if isequaltol(gdata.Bnc         , data.output.Bnc                  , tol) ; disp('ok: gdata.Bnc        ') ; else ; disp('ERROR: gdata.Bnc        '); matching=0; end
if isequaltol(gdata.Bns         , data.output.Bns                  , tol) ; disp('ok: gdata.Bns        ') ; else ; disp('ERROR: gdata.Bns        '); matching=0; end
if isequaltol(gdata.Btemn       , data.output.Btemn                , tol) ; disp('ok: gdata.Btemn      ') ; else ; disp('ERROR: gdata.Btemn      '); matching=0; end
if isequaltol(gdata.Btomn       , data.output.Btomn                , tol) ; disp('ok: gdata.Btomn      ') ; else ; disp('ERROR: gdata.Btomn      '); matching=0; end
if isequaltol(gdata.Bzemn       , data.output.Bzemn                , tol) ; disp('ok: gdata.Bzemn      ') ; else ; disp('ERROR: gdata.Bzemn      '); matching=0; end
if isequaltol(gdata.Bzomn       , data.output.Bzomn                , tol) ; disp('ok: gdata.Bzomn      ') ; else ; disp('ERROR: gdata.Bzomn      '); matching=0; end
if isequaltol(gdata.ForceErr    , data.output.ForceErr             , tol) ; disp('ok: gdata.ForceErr   ') ; else ; disp('ERROR: gdata.ForceErr   '); matching=0; end
if isequaltol(gdata.Igeometry   , data.input.physics.Igeometry     , tol) ; disp('ok: gdata.Igeometry  ') ; else ; disp('ERROR: gdata.Igeometry  '); matching=0; end
if isequaltol(gdata.Istellsym   , data.input.physics.Istellsym     , tol) ; disp('ok: gdata.Istellsym  ') ; else ; disp('ERROR: gdata.Istellsym  '); matching=0; end
if isequaltol(gdata.Ladiabatic  , data.input.physics.Ladiabatic    , tol) ; disp('ok: gdata.Ladiabatic ') ; else ; disp('ERROR: gdata.Ladiabatic '); matching=0; end
if isequaltol(gdata.Lconstraint , data.input.physics.Lconstraint   , tol) ; disp('ok: gdata.Lconstraint') ; else ; disp('ERROR: gdata.Lconstraint'); matching=0; end
if isequaltol(gdata.Lfreebound  , data.input.physics.Lfreebound    , tol) ; disp('ok: gdata.Lfreebound ') ; else ; disp('ERROR: gdata.Lfreebound '); matching=0; end
if isequaltol(gdata.Lperturbed  , data.input.diagnostics.Lperturbed, tol) ; disp('ok: gdata.Lperturbed ') ; else ; disp('ERROR: gdata.Lperturbed '); matching=0; end
if isequaltol(gdata.Lrad        , data.input.physics.Lrad          , tol) ; disp('ok: gdata.Lrad       ') ; else ; disp('ERROR: gdata.Lrad       '); matching=0; end
if isequaltol(gdata.Mpol        , data.input.physics.Mpol          , tol) ; disp('ok: gdata.Mpol       ') ; else ; disp('ERROR: gdata.Mpol       '); matching=0; end
if isequaltol(gdata.Mrad        , data.output.Mrad                 , tol) ; disp('ok: gdata.Mrad       ') ; else ; disp('ERROR: gdata.Mrad       '); matching=0; end
if isequaltol(gdata.Mvol        , data.output.Mvol                 , tol) ; disp('ok: gdata.Mvol       ') ; else ; disp('ERROR: gdata.Mvol       '); matching=0; end
if isequaltol(gdata.Nfp         , data.input.physics.Nfp           , tol) ; disp('ok: gdata.Nfp        ') ; else ; disp('ERROR: gdata.Nfp        '); matching=0; end
if isequaltol(gdata.Ntor        , data.input.physics.Ntor          , tol) ; disp('ok: gdata.Ntor       ') ; else ; disp('ERROR: gdata.Ntor       '); matching=0; end
if isequaltol(gdata.Nvol        , data.input.physics.Nvol          , tol) ; disp('ok: gdata.Nvol       ') ; else ; disp('ERROR: gdata.Nvol       '); matching=0; end
if isequaltol(gdata.Rbc         , data.output.Rbc                  , tol) ; disp('ok: gdata.Rbc        ') ; else ; disp('ERROR: gdata.Rbc        '); matching=0; end
if isequaltol(gdata.Rbs         , data.output.Rbs                  , tol) ; disp('ok: gdata.Rbs        ') ; else ; disp('ERROR: gdata.Rbs        '); matching=0; end
if isequaltol(gdata.TT          , data.output.TT                   , tol) ; disp('ok: gdata.TT         ') ; else ; disp('ERROR: gdata.TT         '); matching=0; end
if isequaltol(gdata.Vnc         , data.output.Vnc                  , tol) ; disp('ok: gdata.Vnc        ') ; else ; disp('ERROR: gdata.Vnc        '); matching=0; end
if isequaltol(gdata.Vns         , data.output.Vns                  , tol) ; disp('ok: gdata.Vns        ') ; else ; disp('ERROR: gdata.Vns        '); matching=0; end
if isequaltol(gdata.Zbc         , data.output.Zbc                  , tol) ; disp('ok: gdata.Zbc        ') ; else ; disp('ERROR: gdata.Zbc        '); matching=0; end
if isequaltol(gdata.Zbs         , data.output.Zbs                  , tol) ; disp('ok: gdata.Zbs        ') ; else ; disp('ERROR: gdata.Zbs        '); matching=0; end
if isequaltol(gdata.adiabatic   , data.output.adiabatic            , tol) ; disp('ok: gdata.adiabatic  ') ; else ; disp('ERROR: gdata.adiabatic  '); matching=0; end
if isequaltol(gdata.curpol      , data.input.physics.curpol        , tol) ; disp('ok: gdata.curpol     ') ; else ; disp('ERROR: gdata.curpol     '); matching=0; end
if isequaltol(gdata.curtor      , data.input.physics.curtor        , tol) ; disp('ok: gdata.curtor     ') ; else ; disp('ERROR: gdata.curtor     '); matching=0; end
if isequaltol(gdata.dpp         , data.input.diagnostics.dpp       , tol) ; disp('ok: gdata.dpp        ') ; else ; disp('ERROR: gdata.dpp        '); matching=0; end
if isequaltol(gdata.dqq         , data.input.diagnostics.dqq       , tol) ; disp('ok: gdata.dqq        ') ; else ; disp('ERROR: gdata.dqq        '); matching=0; end
if isequaltol(gdata.forcetol    , data.input.global.forcetol       , tol) ; disp('ok: gdata.forcetol   ') ; else ; disp('ERROR: gdata.forcetol   '); matching=0; end
if isequaltol(gdata.gamma       , data.input.physics.gamma         , tol) ; disp('ok: gdata.gamma      ') ; else ; disp('ERROR: gdata.gamma      '); matching=0; end
if isequaltol(gdata.helicity    , data.output.helicity             , tol) ; disp('ok: gdata.helicity   ') ; else ; disp('ERROR: gdata.helicity   '); matching=0; end
if isequaltol(gdata.im          , data.output.im                   , tol) ; disp('ok: gdata.im         ') ; else ; disp('ERROR: gdata.im         '); matching=0; end
if isequaltol(gdata.in          , data.output.in                   , tol) ; disp('ok: gdata.in         ') ; else ; disp('ERROR: gdata.in         '); matching=0; end
if isequaltol(gdata.iota        , data.input.physics.iota          , tol) ; disp('ok: gdata.iota       ') ; else ; disp('ERROR: gdata.iota       '); matching=0; end
if isequaltol(gdata.lmns        , data.output.lmns                 , tol) ; disp('ok: gdata.lmns       ') ; else ; disp('ERROR: gdata.lmns       '); matching=0; end
if isequaltol(gdata.lp          , data.input.physics.lp            , tol) ; disp('ok: gdata.lp         ') ; else ; disp('ERROR: gdata.lp         '); matching=0; end
if isequaltol(gdata.lq          , data.input.physics.lq            , tol) ; disp('ok: gdata.lq         ') ; else ; disp('ERROR: gdata.lq         '); matching=0; end
if isequaltol(gdata.mn          , data.output.mn                   , tol) ; disp('ok: gdata.mn         ') ; else ; disp('ERROR: gdata.mn         '); matching=0; end
if isequaltol(gdata.mu          , data.output.mu                   , tol) ; disp('ok: gdata.mu         ') ; else ; disp('ERROR: gdata.mu         '); matching=0; end
if isequaltol(gdata.mupfits     , data.input.physics.mupfits       , tol) ; disp('ok: gdata.mupfits    ') ; else ; disp('ERROR: gdata.mupfits    '); matching=0; end
if isequaltol(gdata.mupftol     , data.input.physics.mupftol       , tol) ; disp('ok: gdata.mupftol    ') ; else ; disp('ERROR: gdata.mupftol    '); matching=0; end
if isequaltol(gdata.oita        , data.input.physics.oita          , tol) ; disp('ok: gdata.oita       ') ; else ; disp('ERROR: gdata.oita       '); matching=0; end
if isequaltol(gdata.pflux       , data.output.pflux                , tol) ; disp('ok: gdata.pflux      ') ; else ; disp('ERROR: gdata.pflux      '); matching=0; end
if isequaltol(gdata.phiedge     , data.input.physics.phiedge       , tol) ; disp('ok: gdata.phiedge    ') ; else ; disp('ERROR: gdata.phiedge    '); matching=0; end
if isequaltol(gdata.pl          , data.input.physics.pl            , tol) ; disp('ok: gdata.pl         ') ; else ; disp('ERROR: gdata.pl         '); matching=0; end
if isequaltol(gdata.pr          , data.input.physics.pr            , tol) ; disp('ok: gdata.pr         ') ; else ; disp('ERROR: gdata.pr         '); matching=0; end
if isequaltol(gdata.pressure    , data.input.physics.pressure      , tol) ; disp('ok: gdata.pressure   ') ; else ; disp('ERROR: gdata.pressure   '); matching=0; end
if isequaltol(gdata.pscale      , data.input.physics.pscale        , tol) ; disp('ok: gdata.pscale     ') ; else ; disp('ERROR: gdata.pscale     '); matching=0; end
if isequaltol(gdata.ql          , data.input.physics.ql            , tol) ; disp('ok: gdata.ql         ') ; else ; disp('ERROR: gdata.ql         '); matching=0; end
if isequaltol(gdata.qr          , data.input.physics.qr            , tol) ; disp('ok: gdata.qr         ') ; else ; disp('ERROR: gdata.qr         '); matching=0; end
if isequaltol(gdata.rp          , data.input.physics.rp            , tol) ; disp('ok: gdata.rp         ') ; else ; disp('ERROR: gdata.rp         '); matching=0; end
if isequaltol(gdata.rq          , data.input.physics.rq            , tol) ; disp('ok: gdata.rq         ') ; else ; disp('ERROR: gdata.rq         '); matching=0; end
if isequaltol(gdata.tflux       , data.output.tflux                , tol) ; disp('ok: gdata.tflux      ') ; else ; disp('ERROR: gdata.tflux      '); matching=0; end
if isequaltol(gdata.volume      , data.output.volume               , tol) ; disp('ok: gdata.volume     ') ; else ; disp('ERROR: gdata.volume     '); matching=0; end

if isequaltol(gdata.Nt          , data.grid.Nt                     , tol) ; disp('ok: gdata.Nt         ') ; else ; disp('ERROR: gdata.Nt         '); matching=0; end
if isequaltol(gdata.Nz          , data.grid.Nz                     , tol) ; disp('ok: gdata.Nz         ') ; else ; disp('ERROR: gdata.Nz         '); matching=0; end
if isequaltol(gdata.Ntz         , data.grid.Ntz                    , tol) ; disp('ok: gdata.Ntz        ') ; else ; disp('ERROR: gdata.Ntz        '); matching=0; end
if isequaltol(gdata.Rij         , data.grid.Rij                    , tol) ; disp('ok: gdata.Rij        ') ; else ; disp('ERROR: gdata.Rij        '); matching=0; end
if isequaltol(gdata.Zij         , data.grid.Zij                    , tol) ; disp('ok: gdata.Zij        ') ; else ; disp('ERROR: gdata.Zij        '); matching=0; end
if isequaltol(gdata.sg          , data.grid.sg                     , tol) ; disp('ok: gdata.sg         ') ; else ; disp('ERROR: gdata.sg         '); matching=0; end
if isequaltol(gdata.BR          , data.grid.BR                     , tol) ; disp('ok: gdata.BR         ') ; else ; disp('ERROR: gdata.BR         '); matching=0; end
if isequaltol(gdata.Bp          , data.grid.Bp                     , tol) ; disp('ok: gdata.Bp         ') ; else ; disp('ERROR: gdata.Bp         '); matching=0; end
if isequaltol(gdata.BZ          , data.grid.BZ                     , tol) ; disp('ok: gdata.BZ         ') ; else ; disp('ERROR: gdata.BZ         '); matching=0; end

% data on rotational transform from .ext.sp.t.0001.dat and so on
if isequaltol(idata.Bnc         , data.output.Bnc                  , tol) ; disp('ok: idata.Bnc        ') ; else ; disp('ERROR: idata.Bnc        '); matching=0; end
if isequaltol(idata.Bns         , data.output.Bns                  , tol) ; disp('ok: idata.Bns        ') ; else ; disp('ERROR: idata.Bns        '); matching=0; end
if isequaltol(idata.Btemn       , data.output.Btemn                , tol) ; disp('ok: idata.Btemn      ') ; else ; disp('ERROR: idata.Btemn      '); matching=0; end
if isequaltol(idata.Btomn       , data.output.Btomn                , tol) ; disp('ok: idata.Btomn      ') ; else ; disp('ERROR: idata.Btomn      '); matching=0; end
if isequaltol(idata.Bzemn       , data.output.Bzemn                , tol) ; disp('ok: idata.Bzemn      ') ; else ; disp('ERROR: idata.Bzemn      '); matching=0; end
if isequaltol(idata.Bzomn       , data.output.Bzomn                , tol) ; disp('ok: idata.Bzomn      ') ; else ; disp('ERROR: idata.Bzomn      '); matching=0; end
if isequaltol(idata.ForceErr    , data.output.ForceErr             , tol) ; disp('ok: idata.ForceErr   ') ; else ; disp('ERROR: idata.ForceErr   '); matching=0; end
if isequaltol(idata.Igeometry   , data.input.physics.Igeometry     , tol) ; disp('ok: idata.Igeometry  ') ; else ; disp('ERROR: idata.Igeometry  '); matching=0; end
if isequaltol(idata.Istellsym   , data.input.physics.Istellsym     , tol) ; disp('ok: idata.Istellsym  ') ; else ; disp('ERROR: idata.Istellsym  '); matching=0; end
if isequaltol(idata.Ladiabatic  , data.input.physics.Ladiabatic    , tol) ; disp('ok: idata.Ladiabatic ') ; else ; disp('ERROR: idata.Ladiabatic '); matching=0; end
if isequaltol(idata.Lconstraint , data.input.physics.Lconstraint   , tol) ; disp('ok: idata.Lconstraint') ; else ; disp('ERROR: idata.Lconstraint'); matching=0; end
if isequaltol(idata.Lfreebound  , data.input.physics.Lfreebound    , tol) ; disp('ok: idata.Lfreebound ') ; else ; disp('ERROR: idata.Lfreebound '); matching=0; end
if isequaltol(idata.Lperturbed  , data.input.diagnostics.Lperturbed, tol) ; disp('ok: idata.Lperturbed ') ; else ; disp('ERROR: idata.Lperturbed '); matching=0; end
if isequaltol(idata.Lrad        , data.input.physics.Lrad          , tol) ; disp('ok: idata.Lrad       ') ; else ; disp('ERROR: idata.Lrad       '); matching=0; end
if isequaltol(idata.Mpol        , data.input.physics.Mpol          , tol) ; disp('ok: idata.Mpol       ') ; else ; disp('ERROR: idata.Mpol       '); matching=0; end
if isequaltol(idata.Mrad        , data.output.Mrad                 , tol) ; disp('ok: idata.Mrad       ') ; else ; disp('ERROR: idata.Mrad       '); matching=0; end
if isequaltol(idata.Mvol        , data.output.Mvol                 , tol) ; disp('ok: idata.Mvol       ') ; else ; disp('ERROR: idata.Mvol       '); matching=0; end
if isequaltol(idata.Nfp         , data.input.physics.Nfp           , tol) ; disp('ok: idata.Nfp        ') ; else ; disp('ERROR: idata.Nfp        '); matching=0; end
if isequaltol(idata.Ntor        , data.input.physics.Ntor          , tol) ; disp('ok: idata.Ntor       ') ; else ; disp('ERROR: idata.Ntor       '); matching=0; end
if isequaltol(idata.Nvol        , data.input.physics.Nvol          , tol) ; disp('ok: idata.Nvol       ') ; else ; disp('ERROR: idata.Nvol       '); matching=0; end
if isequaltol(idata.Rbc         , data.output.Rbc                  , tol) ; disp('ok: idata.Rbc        ') ; else ; disp('ERROR: idata.Rbc        '); matching=0; end
if isequaltol(idata.Rbs         , data.output.Rbs                  , tol) ; disp('ok: idata.Rbs        ') ; else ; disp('ERROR: idata.Rbs        '); matching=0; end
if isequaltol(idata.TT          , data.output.TT                   , tol) ; disp('ok: idata.TT         ') ; else ; disp('ERROR: idata.TT         '); matching=0; end
if isequaltol(idata.Vnc         , data.output.Vnc                  , tol) ; disp('ok: idata.Vnc        ') ; else ; disp('ERROR: idata.Vnc        '); matching=0; end
if isequaltol(idata.Vns         , data.output.Vns                  , tol) ; disp('ok: idata.Vns        ') ; else ; disp('ERROR: idata.Vns        '); matching=0; end
if isequaltol(idata.Zbc         , data.output.Zbc                  , tol) ; disp('ok: idata.Zbc        ') ; else ; disp('ERROR: idata.Zbc        '); matching=0; end
if isequaltol(idata.Zbs         , data.output.Zbs                  , tol) ; disp('ok: idata.Zbs        ') ; else ; disp('ERROR: idata.Zbs        '); matching=0; end
if isequaltol(idata.adiabatic   , data.output.adiabatic            , tol) ; disp('ok: idata.adiabatic  ') ; else ; disp('ERROR: idata.adiabatic  '); matching=0; end
if isequaltol(idata.curpol      , data.input.physics.curpol        , tol) ; disp('ok: idata.curpol     ') ; else ; disp('ERROR: idata.curpol     '); matching=0; end
if isequaltol(idata.curtor      , data.input.physics.curtor        , tol) ; disp('ok: idata.curtor     ') ; else ; disp('ERROR: idata.curtor     '); matching=0; end
if isequaltol(idata.dpp         , data.input.diagnostics.dpp       , tol) ; disp('ok: idata.dpp        ') ; else ; disp('ERROR: idata.dpp        '); matching=0; end
if isequaltol(idata.dqq         , data.input.diagnostics.dqq       , tol) ; disp('ok: idata.dqq        ') ; else ; disp('ERROR: idata.dqq        '); matching=0; end
if isequaltol(idata.forcetol    , data.input.global.forcetol       , tol) ; disp('ok: idata.forcetol   ') ; else ; disp('ERROR: idata.forcetol   '); matching=0; end
if isequaltol(idata.gamma       , data.input.physics.gamma         , tol) ; disp('ok: idata.gamma      ') ; else ; disp('ERROR: idata.gamma      '); matching=0; end
if isequaltol(idata.helicity    , data.output.helicity             , tol) ; disp('ok: idata.helicity   ') ; else ; disp('ERROR: idata.helicity   '); matching=0; end
if isequaltol(idata.im          , data.output.im                   , tol) ; disp('ok: idata.im         ') ; else ; disp('ERROR: idata.im         '); matching=0; end
if isequaltol(idata.in          , data.output.in                   , tol) ; disp('ok: idata.in         ') ; else ; disp('ERROR: idata.in         '); matching=0; end
if isequaltol(idata.lmns        , data.output.lmns                 , tol) ; disp('ok: idata.lmns       ') ; else ; disp('ERROR: idata.lmns       '); matching=0; end
if isequaltol(idata.lp          , data.input.physics.lp            , tol) ; disp('ok: idata.lp         ') ; else ; disp('ERROR: idata.lp         '); matching=0; end
if isequaltol(idata.lq          , data.input.physics.lq            , tol) ; disp('ok: idata.lq         ') ; else ; disp('ERROR: idata.lq         '); matching=0; end
if isequaltol(idata.mn          , data.output.mn                   , tol) ; disp('ok: idata.mn         ') ; else ; disp('ERROR: idata.mn         '); matching=0; end
if isequaltol(idata.mu          , data.output.mu                   , tol) ; disp('ok: idata.mu         ') ; else ; disp('ERROR: idata.mu         '); matching=0; end
if isequaltol(idata.mupfits     , data.input.physics.mupfits       , tol) ; disp('ok: idata.mupfits    ') ; else ; disp('ERROR: idata.mupfits    '); matching=0; end
if isequaltol(idata.mupftol     , data.input.physics.mupftol       , tol) ; disp('ok: idata.mupftol    ') ; else ; disp('ERROR: idata.mupftol    '); matching=0; end
if isequaltol(idata.oita        , data.input.physics.oita          , tol) ; disp('ok: idata.oita       ') ; else ; disp('ERROR: idata.oita       '); matching=0; end
if isequaltol(idata.pflux       , data.output.pflux                , tol) ; disp('ok: idata.pflux      ') ; else ; disp('ERROR: idata.pflux      '); matching=0; end
if isequaltol(idata.phiedge     , data.input.physics.phiedge       , tol) ; disp('ok: idata.phiedge    ') ; else ; disp('ERROR: idata.phiedge    '); matching=0; end
if isequaltol(idata.pl          , data.input.physics.pl            , tol) ; disp('ok: idata.pl         ') ; else ; disp('ERROR: idata.pl         '); matching=0; end
if isequaltol(idata.pr          , data.input.physics.pr            , tol) ; disp('ok: idata.pr         ') ; else ; disp('ERROR: idata.pr         '); matching=0; end
if isequaltol(idata.pressure    , data.input.physics.pressure      , tol) ; disp('ok: idata.pressure   ') ; else ; disp('ERROR: idata.pressure   '); matching=0; end
if isequaltol(idata.pscale      , data.input.physics.pscale        , tol) ; disp('ok: idata.pscale     ') ; else ; disp('ERROR: idata.pscale     '); matching=0; end
if isequaltol(idata.ql          , data.input.physics.ql            , tol) ; disp('ok: idata.ql         ') ; else ; disp('ERROR: idata.ql         '); matching=0; end
if isequaltol(idata.qr          , data.input.physics.qr            , tol) ; disp('ok: idata.qr         ') ; else ; disp('ERROR: idata.qr         '); matching=0; end
if isequaltol(idata.rp          , data.input.physics.rp            , tol) ; disp('ok: idata.rp         ') ; else ; disp('ERROR: idata.rp         '); matching=0; end
if isequaltol(idata.rq          , data.input.physics.rq            , tol) ; disp('ok: idata.rq         ') ; else ; disp('ERROR: idata.rq         '); matching=0; end
if isequaltol(idata.tflux       , data.output.tflux                , tol) ; disp('ok: idata.tflux      ') ; else ; disp('ERROR: idata.tflux      '); matching=0; end
if isequaltol(idata.volume      , data.output.volume               , tol) ; disp('ok: idata.volume     ') ; else ; disp('ERROR: idata.volume     '); matching=0; end

if isequaltol(idata.iota        , data.transform.fiota(:,2)        , tol) ; disp('ok: idata.iota       ') ; else ; disp('ERROR: idata.iota       '); matching=0; end
if isequaltol(idata.sarr        , data.transform.fiota(:,1)'       , tol) ; disp('ok: idata.sarr       ') ; else ; disp('ERROR: idata.sarr       '); matching=0; end

% pdata from field line tracing; from ".ext.sp.P.xxxx.dat" files
if isequaltol(pdata.Bnc         , data.output.Bnc                  , tol) ; disp('ok: pdata.Bnc        ') ; else ; disp('ERROR: pdata.Bnc        '); matching=0; end
if isequaltol(pdata.Bns         , data.output.Bns                  , tol) ; disp('ok: pdata.Bns        ') ; else ; disp('ERROR: pdata.Bns        '); matching=0; end
if isequaltol(pdata.Btemn       , data.output.Btemn                , tol) ; disp('ok: pdata.Btemn      ') ; else ; disp('ERROR: pdata.Btemn      '); matching=0; end
if isequaltol(pdata.Btomn       , data.output.Btomn                , tol) ; disp('ok: pdata.Btomn      ') ; else ; disp('ERROR: pdata.Btomn      '); matching=0; end
if isequaltol(pdata.Bzemn       , data.output.Bzemn                , tol) ; disp('ok: pdata.Bzemn      ') ; else ; disp('ERROR: pdata.Bzemn      '); matching=0; end
if isequaltol(pdata.Bzomn       , data.output.Bzomn                , tol) ; disp('ok: pdata.Bzomn      ') ; else ; disp('ERROR: pdata.Bzomn      '); matching=0; end
if isequaltol(pdata.ForceErr    , data.output.ForceErr             , tol) ; disp('ok: pdata.ForceErr   ') ; else ; disp('ERROR: pdata.ForceErr   '); matching=0; end
if isequaltol(pdata.Igeometry   , data.input.physics.Igeometry     , tol) ; disp('ok: pdata.Igeometry  ') ; else ; disp('ERROR: pdata.Igeometry  '); matching=0; end
if isequaltol(pdata.Istellsym   , data.input.physics.Istellsym     , tol) ; disp('ok: pdata.Istellsym  ') ; else ; disp('ERROR: pdata.Istellsym  '); matching=0; end
if isequaltol(pdata.Ladiabatic  , data.input.physics.Ladiabatic    , tol) ; disp('ok: pdata.Ladiabatic ') ; else ; disp('ERROR: pdata.Ladiabatic '); matching=0; end
if isequaltol(pdata.Lconstraint , data.input.physics.Lconstraint   , tol) ; disp('ok: pdata.Lconstraint') ; else ; disp('ERROR: pdata.Lconstraint'); matching=0; end
if isequaltol(pdata.Lfreebound  , data.input.physics.Lfreebound    , tol) ; disp('ok: pdata.Lfreebound ') ; else ; disp('ERROR: pdata.Lfreebound '); matching=0; end
if isequaltol(pdata.Lperturbed  , data.input.diagnostics.Lperturbed, tol) ; disp('ok: pdata.Lperturbed ') ; else ; disp('ERROR: pdata.Lperturbed '); matching=0; end
if isequaltol(pdata.Lrad        , data.input.physics.Lrad          , tol) ; disp('ok: pdata.Lrad       ') ; else ; disp('ERROR: pdata.Lrad       '); matching=0; end
if isequaltol(pdata.Mpol        , data.input.physics.Mpol          , tol) ; disp('ok: pdata.Mpol       ') ; else ; disp('ERROR: pdata.Mpol       '); matching=0; end
if isequaltol(pdata.Mrad        , data.output.Mrad                 , tol) ; disp('ok: pdata.Mrad       ') ; else ; disp('ERROR: pdata.Mrad       '); matching=0; end
if isequaltol(pdata.Mvol        , data.output.Mvol                 , tol) ; disp('ok: pdata.Mvol       ') ; else ; disp('ERROR: pdata.Mvol       '); matching=0; end
if isequaltol(pdata.Nfp         , data.input.physics.Nfp           , tol) ; disp('ok: pdata.Nfp        ') ; else ; disp('ERROR: pdata.Nfp        '); matching=0; end
if isequaltol(pdata.Ntor        , data.input.physics.Ntor          , tol) ; disp('ok: pdata.Ntor       ') ; else ; disp('ERROR: pdata.Ntor       '); matching=0; end
if isequaltol(pdata.Nvol        , data.input.physics.Nvol          , tol) ; disp('ok: pdata.Nvol       ') ; else ; disp('ERROR: pdata.Nvol       '); matching=0; end
if isequaltol(pdata.Rbc         , data.output.Rbc                  , tol) ; disp('ok: pdata.Rbc        ') ; else ; disp('ERROR: pdata.Rbc        '); matching=0; end
if isequaltol(pdata.Rbs         , data.output.Rbs                  , tol) ; disp('ok: pdata.Rbs        ') ; else ; disp('ERROR: pdata.Rbs        '); matching=0; end
if isequaltol(pdata.TT          , data.output.TT                   , tol) ; disp('ok: pdata.TT         ') ; else ; disp('ERROR: pdata.TT         '); matching=0; end
if isequaltol(pdata.Vnc         , data.output.Vnc                  , tol) ; disp('ok: pdata.Vnc        ') ; else ; disp('ERROR: pdata.Vnc        '); matching=0; end
if isequaltol(pdata.Vns         , data.output.Vns                  , tol) ; disp('ok: pdata.Vns        ') ; else ; disp('ERROR: pdata.Vns        '); matching=0; end
if isequaltol(pdata.Zbc         , data.output.Zbc                  , tol) ; disp('ok: pdata.Zbc        ') ; else ; disp('ERROR: pdata.Zbc        '); matching=0; end
if isequaltol(pdata.Zbs         , data.output.Zbs                  , tol) ; disp('ok: pdata.Zbs        ') ; else ; disp('ERROR: pdata.Zbs        '); matching=0; end
if isequaltol(pdata.adiabatic   , data.output.adiabatic            , tol) ; disp('ok: pdata.adiabatic  ') ; else ; disp('ERROR: pdata.adiabatic  '); matching=0; end
if isequaltol(pdata.curpol      , data.input.physics.curpol        , tol) ; disp('ok: pdata.curpol     ') ; else ; disp('ERROR: pdata.curpol     '); matching=0; end
if isequaltol(pdata.curtor      , data.input.physics.curtor        , tol) ; disp('ok: pdata.curtor     ') ; else ; disp('ERROR: pdata.curtor     '); matching=0; end
if isequaltol(pdata.dpp         , data.input.diagnostics.dpp       , tol) ; disp('ok: pdata.dpp        ') ; else ; disp('ERROR: pdata.dpp        '); matching=0; end
if isequaltol(pdata.dqq         , data.input.diagnostics.dqq       , tol) ; disp('ok: pdata.dqq        ') ; else ; disp('ERROR: pdata.dqq        '); matching=0; end
if isequaltol(pdata.forcetol    , data.input.global.forcetol       , tol) ; disp('ok: pdata.forcetol   ') ; else ; disp('ERROR: pdata.forcetol   '); matching=0; end
if isequaltol(pdata.gamma       , data.input.physics.gamma         , tol) ; disp('ok: pdata.gamma      ') ; else ; disp('ERROR: pdata.gamma      '); matching=0; end
if isequaltol(pdata.helicity    , data.output.helicity             , tol) ; disp('ok: pdata.helicity   ') ; else ; disp('ERROR: pdata.helicity   '); matching=0; end
if isequaltol(pdata.im          , data.output.im                   , tol) ; disp('ok: pdata.im         ') ; else ; disp('ERROR: pdata.im         '); matching=0; end
if isequaltol(pdata.in          , data.output.in                   , tol) ; disp('ok: pdata.in         ') ; else ; disp('ERROR: pdata.in         '); matching=0; end
if isequaltol(pdata.iota        , data.input.physics.iota          , tol) ; disp('ok: pdata.iota       ') ; else ; disp('ERROR: pdata.iota       '); matching=0; end
if isequaltol(pdata.lmns        , data.output.lmns                 , tol) ; disp('ok: pdata.lmns       ') ; else ; disp('ERROR: pdata.lmns       '); matching=0; end
if isequaltol(pdata.lp          , data.input.physics.lp            , tol) ; disp('ok: pdata.lp         ') ; else ; disp('ERROR: pdata.lp         '); matching=0; end
if isequaltol(pdata.lq          , data.input.physics.lq            , tol) ; disp('ok: pdata.lq         ') ; else ; disp('ERROR: pdata.lq         '); matching=0; end
if isequaltol(pdata.mn          , data.output.mn                   , tol) ; disp('ok: pdata.mn         ') ; else ; disp('ERROR: pdata.mn         '); matching=0; end
if isequaltol(pdata.mu          , data.output.mu                   , tol) ; disp('ok: pdata.mu         ') ; else ; disp('ERROR: pdata.mu         '); matching=0; end
if isequaltol(pdata.mupfits     , data.input.physics.mupfits       , tol) ; disp('ok: pdata.mupfits    ') ; else ; disp('ERROR: pdata.mupfits    '); matching=0; end
if isequaltol(pdata.mupftol     , data.input.physics.mupftol       , tol) ; disp('ok: pdata.mupftol    ') ; else ; disp('ERROR: pdata.mupftol    '); matching=0; end
if isequaltol(pdata.oita        , data.input.physics.oita          , tol) ; disp('ok: pdata.oita       ') ; else ; disp('ERROR: pdata.oita       '); matching=0; end
if isequaltol(pdata.pflux       , data.output.pflux                , tol) ; disp('ok: pdata.pflux      ') ; else ; disp('ERROR: pdata.pflux      '); matching=0; end
if isequaltol(pdata.phiedge     , data.input.physics.phiedge       , tol) ; disp('ok: pdata.phiedge    ') ; else ; disp('ERROR: pdata.phiedge    '); matching=0; end
if isequaltol(pdata.pl          , data.input.physics.pl            , tol) ; disp('ok: pdata.pl         ') ; else ; disp('ERROR: pdata.pl         '); matching=0; end
if isequaltol(pdata.pr          , data.input.physics.pr            , tol) ; disp('ok: pdata.pr         ') ; else ; disp('ERROR: pdata.pr         '); matching=0; end
if isequaltol(pdata.pressure    , data.input.physics.pressure      , tol) ; disp('ok: pdata.pressure   ') ; else ; disp('ERROR: pdata.pressure   '); matching=0; end
if isequaltol(pdata.pscale      , data.input.physics.pscale        , tol) ; disp('ok: pdata.pscale     ') ; else ; disp('ERROR: pdata.pscale     '); matching=0; end
if isequaltol(pdata.ql          , data.input.physics.ql            , tol) ; disp('ok: pdata.ql         ') ; else ; disp('ERROR: pdata.ql         '); matching=0; end
if isequaltol(pdata.qr          , data.input.physics.qr            , tol) ; disp('ok: pdata.qr         ') ; else ; disp('ERROR: pdata.qr         '); matching=0; end
if isequaltol(pdata.rp          , data.input.physics.rp            , tol) ; disp('ok: pdata.rp         ') ; else ; disp('ERROR: pdata.rp         '); matching=0; end
if isequaltol(pdata.rq          , data.input.physics.rq            , tol) ; disp('ok: pdata.rq         ') ; else ; disp('ERROR: pdata.rq         '); matching=0; end
if isequaltol(pdata.tflux       , data.output.tflux                , tol) ; disp('ok: pdata.tflux      ') ; else ; disp('ERROR: pdata.tflux      '); matching=0; end
if isequaltol(pdata.volume      , data.output.volume               , tol) ; disp('ok: pdata.volume     ') ; else ; disp('ERROR: pdata.volume     '); matching=0; end

if isequaltol(pdata.npoinc      , data.input.diagnostics.nPpts     , tol) ; disp('ok: pdata.npoinc     ') ; else ; disp('ERROR: pdata.npoinc     '); matching=0; end
if isequaltol(pdata.th_lines    , data.poincare.t                  , tol) ; disp('ok: pdata.th_lines   ') ; else ; disp('ERROR: pdata.th_lines   '); matching=0; end
if isequaltol(pdata.rho_lines   , data.poincare.rho                , tol) ; disp('ok: pdata.rho_lines  ') ; else ; disp('ERROR: pdata.rho_lines  '); matching=0; end
if isequaltol(pdata.R_lines     , data.poincare.R                  , tol) ; disp('ok: pdata.R_lines    ') ; else ; disp('ERROR: pdata.R_lines    '); matching=0; end
if isequaltol(pdata.Z_lines     , data.poincare.Z                  , tol) ; disp('ok: pdata.npoinc     ') ; else ; disp('ERROR: pdata.Z_lines    '); matching=0; end

if (matching == 0)
  disp('Not maching :(')
else
  disp('Matching :)')
end

end

