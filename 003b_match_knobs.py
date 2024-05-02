from cpymad.madx import Madx
import xtrack as xt
import numpy as np

mad = Madx()

mad.input('''
    call, file="acc-models-lhc/lhc.seq";

    beam, sequence=lhcb1, bv= 1, particle=proton, charge=1, mass=0.938272046,
    pc= 450.0,   npart=1.2e11,kbunch=2556, ex=5.2126224777777785e-09,ey=5.2126224777777785e-09;
    beam, sequence=lhcb2, bv=-1, particle=proton, charge=1, mass=0.938272046,
    pc= 450.0,   npart=1.2e11,kbunch=2556, ex=5.2126224777777785e-09,ey=5.2126224777777785e-09;

    call,file="acc-models-lhc/operation/optics/R2024aRP_A11mC11mA10mL10m_PhaseKnob100ON.madx";
''')

twmad_b1  = mad.twiss(sequence='lhcb1', table='twb1')
twmad_b2  = mad.twiss(sequence='lhcb2', table='twb2')

collider = xt.Multiline.from_madx(madx=mad)
collider.lhcb1.twiss_default['method'] = '4d'
collider.lhcb2.twiss_default['method'] = '4d'
collider.lhcb2.twiss_default['reverse'] = True

# All orbit correctors off
corrector_names = collider.vars.get_table().rows['acb.*'].name
for nn in corrector_names:
    collider.vars[nn] = 0.0

twflat = collider.twiss()

assert_allclose = np.testing.assert_allclose
assert_allclose(twflat.lhcb1.x, 0.0, atol=1e-14)
assert_allclose(twflat.lhcb1.y, 0.0, atol=1e-14)
assert_allclose(twflat.lhcb2.x, 0.0, atol=1e-14)
assert_allclose(twflat.lhcb2.y, 0.0, atol=1e-14)