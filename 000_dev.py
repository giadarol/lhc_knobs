from cpymad.madx import Madx
import xtrack as xt

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

# Match tune knob
b1_mqt_circuits = [
 'kqtf.a12b1',
 'kqtf.a23b1',
 'kqtf.a34b1',
 'kqtf.a45b1',
 'kqtf.a56b1',
 'kqtf.a67b1',
 'kqtf.a78b1',
 'kqtf.a81b1']

b2_mqt_circuits = [
 'kqtf.a12b2',
 'kqtf.a23b2',
 'kqtf.a34b2',
 'kqtf.a45b2',
 'kqtf.a56b2',
 'kqtf.a67b2',
 'kqtf.a78b2',
 'kqtf.a81b2']

for nn in b1_mqt_circuits + b2_mqt_circuits:
    collider.vars['old_' + nn] = collider.vars[nn]._value
    collider.vars[nn] = collider.vars[nn]._value

dq_match = 1e-4
twb1 = collider.lhcb1.twiss()

opt_b1_dqx = collider.lhcb1.match_knob(
    knob_name='dqx.b1',
    knob_value_start=0.,
    knob_value_end=dq_match,
    run=False,
    vary=xt.VaryList(b1_mqt_circuits, step=1e-8),
    targets=[
        twb1.target('qx', twb1.qx + dq_match, tolerance=1e-7),
        twb1.target('qy', twb1.qy, tolerance=1e-7)
    ]
)

opt_b1_dqx.solve()
opt_b1_dqx.generate_knob()



