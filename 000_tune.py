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

mqt_circuits={}
mqt_circuits['b1'] = [
 'kqtf.a12b1',
 'kqtf.a23b1',
 'kqtf.a34b1',
 'kqtf.a45b1',
 'kqtf.a56b1',
 'kqtf.a67b1',
 'kqtf.a78b1',
 'kqtf.a81b1',
 'kqtd.a12b1',
 'kqtd.a23b1',
 'kqtd.a34b1',
 'kqtd.a45b1',
 'kqtd.a56b1',
 'kqtd.a67b1',
 'kqtd.a78b1',
 'kqtd.a81b1']

mqt_circuits['b2'] = [
 'kqtf.a12b2',
 'kqtf.a23b2',
 'kqtf.a34b2',
 'kqtf.a45b2',
 'kqtf.a56b2',
 'kqtf.a67b2',
 'kqtf.a78b2',
 'kqtf.a81b2',
 'kqtd.a12b2',
 'kqtd.a23b2',
 'kqtd.a34b2',
 'kqtd.a45b2',
 'kqtd.a56b2',
 'kqtd.a67b2',
 'kqtd.a78b2',
 'kqtd.a81b2']

for nn in mqt_circuits['b1'] + mqt_circuits['b2']:
    collider.vars['old_' + nn] = collider.vars[nn]._expr
    collider.vars[nn] = collider.vars[nn]._value

optimizers = {'b1': {}, 'b2': {}}
dq_match = 1e-3
for bname in ['b1', 'b2']:
    tw = collider[f'lhc{bname}'].twiss()
    opt_qx = collider[f'lhc{bname}'].match_knob(
        knob_name=f'dqx.{bname}',
        knob_value_start=0.,
        knob_value_end=dq_match,
        run=False,
        vary=xt.VaryList(mqt_circuits[bname], step=1e-8),
        targets=[
            tw.target('qx', tw.qx + dq_match, tol=1e-7),
            tw.target('qy', tw.qy, tol=1e-7)
        ]
    )

    opt_qx.solve()
    opt_qx.generate_knob()
    optimizers[bname]['qx'] = opt_qx

    opt_qy = collider[f'lhc{bname}'].match_knob(
        knob_name=f'dqy.{bname}',
        knob_value_start=0.,
        knob_value_end=dq_match,
        run=False,
        vary=xt.VaryList(mqt_circuits[bname], step=1e-8),
        targets=[
            tw.target('qx', tw.qx, tolerance=1e-7),
            tw.target('qy', tw.qy + dq_match, tolerance=1e-7)
        ]
    )

    opt_qy.solve()
    opt_qy.generate_knob()
    optimizers[bname]['qy'] = opt_qy
