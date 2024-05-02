from cpymad.madx import Madx
import xtrack as xt
import numpy as np
import json

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

# Load configs
with open('ip_orbit_knobs_configs.json', 'r') as fid:
    configs = json.load(fid)

nn = list(configs.keys())[0] # To be replaced by a loop over all knobs

conf = configs[nn]
ipn = conf['ip']

# Build targets
targets = []
for lname in ['lhcb1', 'lhcb2']:
    for tname in conf['targets'][lname]:
        targets.append(xt.Target(tname, line=lname, at='ip'+str(ipn),
                                 value=conf['targets'][lname][tname]))

start_b1 = f'e.ds.l{ipn}.b1'
start_b2 = f'e.ds.l{ipn}.b2'
end_b1 = f's.ds.r{ipn}.b1'
end_b2 = f's.ds.r{ipn}.b2'

plane = conf['plane']

targets += [
    xt.Target(plane, 0, line='lhcb1', at=xt.END),
    xt.Target(plane, 0, line='lhcb2', at=xt.END),
    xt.Target('p' + plane, 0, line='lhcb2', at=xt.END),
    xt.Target('p' + plane, 0, line='lhcb1', at=xt.END),
]

vary = []
for cc in conf['correctors']:
    vary.append(xt.Vary(cc, step=1e-8))

