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

# Extract list of circuits controlled by cromaticity knobs


ms_circuits = {}
ms_circuits['b1'] = [
 'ksf1.a12b1',
 'ksf1.a23b1',
 'ksf1.a34b1',
 'ksf1.a45b1',
 'ksf1.a56b1',
 'ksf1.a67b1',
 'ksf1.a78b1',
 'ksf1.a81b1',
 'ksf2.a12b1',
 'ksf2.a23b1',
 'ksf2.a34b1',
 'ksf2.a45b1',
 'ksf2.a56b1',
 'ksf2.a67b1',
 'ksf2.a78b1',
 'ksf2.a81b1',
 'ksd1.a12b1',
 'ksd1.a23b1',
 'ksd1.a34b1',
 'ksd1.a45b1',
 'ksd1.a56b1',
 'ksd1.a67b1',
 'ksd1.a78b1',
 'ksd1.a81b1',
 'ksd2.a12b1',
 'ksd2.a23b1',
 'ksd2.a34b1',
 'ksd2.a45b1',
 'ksd2.a56b1',
 'ksd2.a67b1',
 'ksd2.a78b1',
 'ksd2.a81b1']

ms_circuits['b2'] = [
 'ksf1.a12b2',
 'ksf1.a23b2',
 'ksf1.a34b2',
 'ksf1.a45b2',
 'ksf1.a56b2',
 'ksf1.a67b2',
 'ksf1.a78b2',
 'ksf1.a81b2',
 'ksf2.a12b2',
 'ksf2.a23b2',
 'ksf2.a34b2',
 'ksf2.a45b2',
 'ksf2.a56b2',
 'ksf2.a67b2',
 'ksf2.a78b2',
 'ksf2.a81b2',
 'ksd1.a12b2',
 'ksd1.a23b2',
 'ksd1.a34b2',
 'ksd1.a45b2',
 'ksd1.a56b2',
 'ksd1.a67b2',
 'ksd1.a78b2',
 'ksd1.a81b2',
 'ksd2.a12b2',
 'ksd2.a23b2',
 'ksd2.a34b2',
 'ksd2.a45b2',
 'ksd2.a56b2',
 'ksd2.a67b2',
 'ksd2.a78b2',
 'ksd2.a81b2']


for nn in ms_circuits['b1'] + ms_circuits['b2']:
    collider.vars['old_' + nn] = collider.vars[nn]._expr
    collider.vars[nn] = collider.vars[nn]._value

optimizers = {'b1': {}, 'b2': {}}
d_chrom_match = 0.5
for bname in ['b1', 'b2']:
    tw = collider[f'lhc{bname}'].twiss()
    opt_qpx = collider[f'lhc{bname}'].match_knob(
        knob_name=f'dqpx.{bname}_op',
        knob_value_start=0.,
        knob_value_end=d_chrom_match,
        run=False,
        vary=xt.VaryList(ms_circuits[bname], step=1e-5),
        targets=[
            xt.Target('dqx', tw.dqx + d_chrom_match, tol=1e-4),
            xt.Target('dqy', tw.dqy, tol=1e-4)
        ]
    )

    opt_qpx.solve()
    opt_qpx.generate_knob()
    optimizers[bname]['qpy'] = opt_qpx

    opt_qpy = collider[f'lhc{bname}'].match_knob(
        knob_name=f'dqpy.{bname}_op',
        knob_value_start=0.,
        knob_value_end=d_chrom_match,
        run=False,
        vary=xt.VaryList(ms_circuits[bname], step=1e-5),
        targets=[
            xt.Target('dqx', tw.dqx, tol=1e-4),
            xt.Target('dqy', tw.dqy + d_chrom_match, tol=1e-4)
        ]
    )

    opt_qpy.solve()
    opt_qpy.generate_knob()
    optimizers[bname]['qpy'] = opt_qpy

# Test the knobs

# Correct to zero
collider.lhcb1.match(
    vary=xt.VaryList(['dqpx.b1_op', 'dqpy.b1_op'], step=1e-4),
    targets=xt.TargetSet(dqx=0., dqy=0., tol=1e-4))
collider.lhcb2.match(
    vary=xt.VaryList(['dqpx.b2_op', 'dqpy.b2_op'], step=1e-4),
    targets=xt.TargetSet(dqx=0., dqy=0., tol=1e-4))

# Apply deltas
collider.vars['dqpx.b1_op'] += 2
collider.vars['dqpy.b1_op'] += 4
collider.vars['dqpx.b2_op'] += 3
collider.vars['dqpy.b2_op'] += 5

twtest = collider.twiss()