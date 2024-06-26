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

print('cmin off: ')
print('dqmin b1:', twmad_b1.summary['dqmin'])
print('dqmin b2:', twmad_b2.summary['dqmin'])

mad.globals['cmrs.b1_op'] = 1e-3
mad.globals['cmrs.b2_op'] = 1e-3

twmad_b1_cm  = mad.twiss(sequence='lhcb1', table='twb1')
twmad_b2_cm  = mad.twiss(sequence='lhcb2', table='twb2')

print('cmin on: ')
print('dqmin b1:', twmad_b1_cm.summary['dqmin'])
print('dqmin b2:', twmad_b2_cm.summary['dqmin'])