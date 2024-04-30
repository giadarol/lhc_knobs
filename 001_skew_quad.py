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

# Match tune knob
mqs_circuits_4_quads={}
mqs_circuits_2_quads={}
mqs_circuits_4_quads['b1'] = [
 'kqs.a23b1',
 'kqs.a45b1',
 'kqs.a67b1',
 'kqs.a81b1',
]
mqs_circuits_2_quads['b1'] = [
 'kqs.l2b1',
 'kqs.l4b1',
 'kqs.l6b1',
 'kqs.l8b1',
 'kqs.r1b1',
 'kqs.r3b1',
 'kqs.r5b1',
 'kqs.r7b1']

mqs_circuits_4_quads['b2'] = [
 'kqs.a12b2',
 'kqs.a34b2',
 'kqs.a56b2',
 'kqs.a78b2']

mqs_circuits_2_quads['b2'] = [
 'kqs.l1b2',
 'kqs.l3b2',
 'kqs.l5b2',
 'kqs.l7b2',
 'kqs.r2b2',
 'kqs.r4b2',
 'kqs.r6b2',
 'kqs.r8b2']

# see Eq. 47 in https://cds.cern.ch/record/522049/files/lhc-project-report-501.pdf
class ActionCmin(xt.Action):
    def __init__(self, line):
        self.line = line
    def run(self):
        tw = self.line.twiss(strengths=True)
        k1sl = tw['k1sl']
        c_min = 1 / (2*np.pi) * np.sum(k1sl * np.sqrt(tw.betx * tw.bety)
                                * np.exp(1j * 2 * np.pi * (tw.mux - tw.muy)))
        return {'c_min_re': c_min.real, 'c_min_im': c_min.imag}

act_cmin_b1 = ActionCmin(collider['lhcb1'])
act_cmin_b2 = ActionCmin(collider['lhcb2'])

for nn in (mqs_circuits_4_quads['b1'] + mqs_circuits_2_quads['b1']
           + mqs_circuits_4_quads['b2'] + mqs_circuits_2_quads['b2']):
    collider.vars['old_' + nn] = collider.vars[nn]._value
    collider.vars[nn] = collider.vars[nn]._value

optimizers = {'b1': {}, 'b2': {}}

c_min_match = 1e-4
for bname in ['b1', 'b2']:
    line = collider[f'lhc{bname}']
    act_cmin = ActionCmin(line)

    assert np.abs(act_cmin.run()['c_min_re']) < 1e-6
    assert np.abs(act_cmin.run()['c_min_im']) < 1e-6

    opt_re = line.match_knob(
        run=False,
        knob_name=f'cmrs.{bname}_op',
        knob_value_start=0,
        knob_value_end=c_min_match,
        vary=[xt.VaryList(mqs_circuits_2_quads[bname], step=5e-5),
            xt.VaryList(mqs_circuits_4_quads[bname], step=5e-5, weight=2)],
        targets=[
            act_cmin.target('c_min_re', value=c_min_match, tol=1e-8),
            act_cmin.target('c_min_im', value=0, tol=1e-8),
        ])
    opt_re.solve()
    opt_re.generate_knob()
    optimizers[bname]['re'] = opt_re

    opt_im = line.match_knob(
        run=False,
        knob_name=f'cmis.{bname}_op',
        knob_value_start=0,
        knob_value_end=c_min_match,
        vary=[xt.VaryList(mqs_circuits_2_quads[bname], step=5e-5),
            xt.VaryList(mqs_circuits_4_quads[bname], step=5e-5, weight=2)],
        targets=[
            act_cmin.target('c_min_re', value=0, tol=1e-8),
            act_cmin.target('c_min_im', value=c_min_match, tol=1e-8),
        ])
    opt_im.solve()
    opt_im.generate_knob()
    optimizers[bname]['im'] = opt_im

