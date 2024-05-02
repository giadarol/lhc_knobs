import xtrack as xt

dct = {'on_sep1_h': {'acbch5.l1b2': -1.41918166577e-05,
  'acbch5.r1b1': 1.64239081776e-05,
  'acbch6.l1b1': -3.4974463131e-06,
  'acbch6.r1b2': 3.58269623014e-06,
  'acbxh1.l1': 8e-06,
  'acbxh1.r1': 8e-06,
  'acbxh2.l1': 8e-06,
  'acbxh2.r1': 8e-06,
  'acbxh3.l1': 8e-06,
  'acbxh3.r1': 8e-06,
  'acbyhs4.l1b1': 1.36297113998e-05,
  'acbyhs4.l1b2': 9.3708594619e-06,
  'acbyhs4.r1b1': -1.24156012491e-05,
  'acbyhs4.r1b2': -1.39157673175e-05},
 'on_sep2h': {'acbch5.r2b2': 0.0,
  'acbchs5.r2b1': 7.205639217199999e-06,
  'acbchs5.r2b2': 5.39723817366e-06,
  'acbxh1.l2': 0.0,
  'acbxh1.r2': 9e-06,
  'acbxh2.l2': 1.35e-05,
  'acbxh2.r2': 9e-06,
  'acbxh3.l2': 1.35e-05,
  'acbxh3.r2': 9e-06,
  'acbyh4.l2b2': 0.0,
  'acbyh4.r2b1': 0.0,
  'acbyh5.l2b1': 0.0,
  'acbyhs4.l2b1': 8.03040506369e-06,
  'acbyhs4.l2b2': -6.55793522704e-06,
  'acbyhs4.r2b1': 6.626871073399999e-06,
  'acbyhs4.r2b2': -1.8905667239800002e-05,
  'acbyhs5.l2b1': -1.5501793291099998e-06,
  'acbyhs5.l2b2': -6.9773770495e-06},
 'on_sep5_v': {'acbcv5.l5b1': -1.02875467053e-05,
  'acbcv5.r5b2': 9.06319153047e-06,
  'acbcv6.l5b2': -1.05495498088e-05,
  'acbcv6.r5b1': 9.86349425437e-06,
  'acbxv1.l5': 8e-06,
  'acbxv1.r5': 8e-06,
  'acbxv2.l5': 8e-06,
  'acbxv2.r5': 8e-06,
  'acbxv3.l5': 8e-06,
  'acbxv3.r5': 8e-06,
  'acbyvs4.l5b1': 2.32108182498e-05,
  'acbyvs4.l5b2': 4.25628214283e-06,
  'acbyvs4.r5b1': -2.39208859279e-06,
  'acbyvs4.r5b2': -2.10127379873e-05},
 'on_sep8v': {'acbcv5.l8b2': 0.0,
  'acbcvs5.l8b1': 4.704301035009999e-06,
  'acbcvs5.l8b2': 4.59473378142e-06,
  'acbxv1.l8': 9e-06,
  'acbxv1.r8': 9e-06,
  'acbxv2.l8': 9e-06,
  'acbxv2.r8': 9e-06,
  'acbxv3.l8': 9e-06,
  'acbxv3.r8': 9e-06,
  'acbyv4.l8b1': 0.0,
  'acbyv4.r8b2': 0.0,
  'acbyv5.r8b1': 0.0,
  'acbyvs4.l8b1': 8.407397301870002e-06,
  'acbyvs4.l8b2': -1.69056856716e-05,
  'acbyvs4.r8b1': 1.745385752e-05,
  'acbyvs4.r8b2': -8.443972688370001e-06,
  'acbyvs5.r8b1': -4.33151486387e-06,
  'acbyvs5.r8b2': -4.72205870866e-06},
 'on_sep1_v': {'acbcv5.l1b1': -1.41681624872e-05,
  'acbcv5.r1b2': 1.64237024066e-05,
  'acbcv6.l1b2': -3.47632625045e-06,
  'acbcv6.r1b1': 3.6084306567e-06,
  'acbxv1.l1': -8e-06,
  'acbxv1.r1': -8e-06,
  'acbxv2.l1': -8e-06,
  'acbxv2.r1': -8e-06,
  'acbxv3.l1': -8e-06,
  'acbxv3.r1': -8e-06,
  'acbyvs4.l1b1': 9.37344313262e-06,
  'acbyvs4.l1b2': 1.36566804063e-05,
  'acbyvs4.r1b1': -1.38925278577e-05,
  'acbyvs4.r1b2': -1.24152779217e-05},
 'on_sep2v': {'acbcv5.r2b1': 0.0,
  'acbcvs5.r2b1': -5.43938795831e-06,
  'acbcvs5.r2b2': -7.28659552868e-06,
  'acbxv1.l2': 9e-06,
  'acbxv1.r2': 9e-06,
  'acbxv2.l2': 9e-06,
  'acbxv2.r2': 9e-06,
  'acbxv3.l2': 9e-06,
  'acbxv3.r2': 9e-06,
  'acbyv4.l2b1': 0.0,
  'acbyv4.r2b2': 0.0,
  'acbyv5.l2b2': 0.0,
  'acbyvs4.l2b1': 9.67871365304e-06,
  'acbyvs4.l2b2': -1.69729431706e-05,
  'acbyvs4.r2b1': 1.87114720092e-05,
  'acbyvs4.r2b2': -6.93506266163e-06,
  'acbyvs5.l2b1': 4.8484845704e-06,
  'acbyvs5.l2b2': 4.34181842841e-06},
 'on_sep5_h': {'acbch5.l5b2': -1.03447685685e-05,
  'acbch5.r5b1': 9.05992577644e-06,
  'acbch6.l5b1': -1.04146390652e-05,
  'acbch6.r5b2': 1.00099005116e-05,
  'acbxh1.l5': -8e-06,
  'acbxh1.r5': -8e-06,
  'acbxh2.l5': -8e-06,
  'acbxh2.r5': -8e-06,
  'acbxh3.l5': -8e-06,
  'acbxh3.r5': -8e-06,
  'acbyhs4.l5b1': 3.13282618249e-06,
  'acbyhs4.l5b2': 2.28458183833e-05,
  'acbyhs4.r5b1': -2.13925054414e-05,
  'acbyhs4.r5b2': -3.49531196308e-06},
 'on_sep8h': {'acbch5.l8b1': 0.0,
  'acbchs5.l8b1': -4.40800786381e-06,
  'acbchs5.l8b2': -4.57801904082e-06,
  'acbxh1.l8': 9e-06,
  'acbxh1.r8': 9e-06,
  'acbxh2.l8': 9e-06,
  'acbxh2.r8': 9e-06,
  'acbxh3.l8': 9e-06,
  'acbxh3.r8': 9e-06,
  'acbyh4.l8b2': 0.0,
  'acbyh4.r8b1': 0.0,
  'acbyh5.r8b2': 0.0,
  'acbyhs4.l8b1': 1.69035964274e-05,
  'acbyhs4.l8b2': -8.43572526061e-06,
  'acbyhs4.r8b1': 9.30146960875e-06,
  'acbyhs4.r8b2': -1.69260077017e-05,
  'acbyhs5.r8b1': 4.83652230929e-06,
  'acbyhs5.r8b2': 4.38410426426e-06},
 'on_x1_v': {'acbcv5.l1b1': -4.49892171021e-08,
  'acbcv5.r1b2': -5.2151408641e-08,
  'acbcv6.l1b2': 4.75980448766e-08,
  'acbcv6.r1b1': 4.94068274744e-08,
  'acbxv1.l1': 5.49019607843e-08,
  'acbxv1.r1': -5.49019607843e-08,
  'acbxv2.l1': 5.49019607843e-08,
  'acbxv2.r1': -5.49019607843e-08,
  'acbxv3.l1': 5.49019607843e-08,
  'acbxv3.r1': -5.49019607843e-08,
  'acbyvs4.l1b1': -2.17704181511e-07,
  'acbyvs4.l1b2': 2.93325423015e-07,
  'acbyvs4.r1b1': 2.90096184016e-07,
  'acbyvs4.r1b2': -2.08045209558e-07},
 'on_x2v': {'acbcv5.r2b1': 0.0,
  'acbcvs5.r2b1': 2.63924957448e-07,
  'acbcvs5.r2b2': -5.5611394134e-08,
  'acbxv1.l2': 5.88235294118e-09,
  'acbxv1.r2': -5.88235294118e-09,
  'acbxv2.l2': 5.88235294118e-09,
  'acbxv2.r2': -5.88235294118e-09,
  'acbxv3.l2': 5.88235294118e-09,
  'acbxv3.r2': -5.88235294118e-09,
  'acbyv4.l2b1': 0.0,
  'acbyv4.r2b2': 0.0,
  'acbyv5.l2b2': 0.0,
  'acbyvs4.l2b1': -3.59306469082e-07,
  'acbyvs4.l2b2': 5.77672138807e-08,
  'acbyvs4.r2b1': -2.65881050814e-08,
  'acbyvs4.r2b2': -3.38366891365e-07,
  'acbyvs5.l2b1': -3.70036866051e-08,
  'acbyvs5.l2b2': 2.10669686866e-07},
 'on_x5_h': {'acbch5.l5b2': 1.02486278372e-07,
  'acbch5.r5b1': 8.97572884285e-08,
  'acbch6.l5b1': -1.6842005584e-07,
  'acbch6.r5b2': -1.61874731447e-07,
  'acbxh1.l5': 5.49019607843e-08,
  'acbxh1.r5': -5.49019607843e-08,
  'acbxh2.l5': 5.49019607843e-08,
  'acbxh2.r5': -5.49019607843e-08,
  'acbxh3.l5': 5.49019607843e-08,
  'acbxh3.r5': -5.49019607843e-08,
  'acbyhs4.l5b1': -9.83004449266e-08,
  'acbyhs4.l5b2': 2.24549528623e-07,
  'acbyhs4.r5b1': 2.38947550477e-07,
  'acbyhs4.r5b2': -9.24385638824e-08},
 'on_x1_h': {'acbch5.l1b2': 4.506876961269906e-08,
  'acbch5.r1b1': 5.2157133283597154e-08,
  'acbch6.l1b1': -4.7890510938800164e-08,
  'acbch6.r1b2': -4.9057882341900174e-08,
  'acbxh1.l1': 5.490196078429958e-08,
  'acbxh1.r1': -5.490196078429958e-08,
  'acbxh2.l1': 5.490196078429958e-08,
  'acbxh2.r1': -5.490196078429958e-08,
  'acbxh3.l1': 5.490196078429958e-08,
  'acbxh3.r1': -5.490196078429958e-08,
  'acbyhs4.l1b1': -2.936865112280011e-07,
  'acbyhs4.l1b2': 2.177064870590009e-07,
  'acbyhs4.r1b1': 2.0803740038799836e-07,
  'acbyhs4.r1b2': -2.89769479559999e-07},
 'on_x2h': {'acbch5.r2b2': 0.0,
  'acbchs5.r2b1': 5.500712482970298e-08,
  'acbchs5.r2b2': -2.619149321510001e-07,
  'acbxh1.l2': 0.0,
  'acbxh1.r2': -5.882352941177012e-09,
  'acbxh2.l2': 8.82352941176213e-09,
  'acbxh2.r2': -5.882352941177012e-09,
  'acbxh3.l2': 8.82352941176213e-09,
  'acbxh3.r2': -5.882352941177012e-09,
  'acbyh4.l2b2': 0.0,
  'acbyh4.r2b1': 0.0,
  'acbyh5.l2b1': 0.0,
  'acbyhs4.l2b1': -3.016846255819841e-08,
  'acbyhs4.l2b2': 3.516576042850015e-07,
  'acbyhs4.r2b1': 3.360059857369998e-07,
  'acbyhs4.r2b2': 3.609186240289533e-08,
  'acbyhs5.l2b1': -2.0978544129400017e-07,
  'acbyhs5.l2b2': 3.489932478159961e-08},
 'on_x5_v': {'acbcv5.l5b1': -1.0191059306399972e-07,
  'acbcv5.r5b2': -8.978187541590162e-08,
  'acbcv6.l5b2': 1.7060064662799925e-07,
  'acbcv6.r5b1': 1.5950619194600115e-07,
  'acbxv1.l5': 5.490196078429958e-08,
  'acbxv1.r5': -5.490196078429958e-08,
  'acbxv2.l5': 5.490196078429958e-08,
  'acbxv2.r5': -5.490196078429958e-08,
  'acbxv3.l5': 5.490196078429958e-08,
  'acbxv3.r5': -5.490196078429958e-08,
  'acbyvs4.l5b1': -2.209491956449975e-07,
  'acbyvs4.l5b2': 8.013405968279948e-08,
  'acbyvs4.r5b1': 1.1028061799399985e-07,
  'acbyvs4.r5b2': -2.427238375419971e-07},
 'on_x8v': {'acbcv5.l8b2': 0.0,
  'acbcvs5.l8b1': 1.0369522268200402e-08,
  'acbcvs5.l8b2': 1.956592257500014e-07,
  'acbxv1.l8': 1.1764705882401458e-08,
  'acbxv1.r8': -1.1764705882401458e-08,
  'acbxv2.l8': 1.1764705882401458e-08,
  'acbxv2.r8': -1.1764705882401458e-08,
  'acbxv3.l8': 1.1764705882401458e-08,
  'acbxv3.r8': -1.1764705882401458e-08,
  'acbyv4.l8b1': 0.0,
  'acbyv4.r8b2': 0.0,
  'acbyv5.r8b1': 0.0,
  'acbyvs4.l8b1': -3.705986005320014e-07,
  'acbyvs4.l8b2': 1.1718075751199826e-07,
  'acbyvs4.r8b1': 9.383772082330024e-08,
  'acbyvs4.r8b2': -3.705179788330009e-07,
  'acbyvs5.r8b1': 1.8445049648999976e-07,
  'acbyvs5.r8b2': 1.040864846200138e-08},
 'on_a2': {'acbch5.r2b2': 0.0,
  'acbchs5.r2b1': 1.2587850689e-07,
  'acbchs5.r2b2': 2.93963812481e-07,
  'acbyh4.l2b2': 0.0,
  'acbyh4.r2b1': 0.0,
  'acbyh5.l2b1': 0.0,
  'acbyhs4.l2b1': 3.68242425339e-08,
  'acbyhs4.l2b2': -3.30398196602e-07,
  'acbyhs4.r2b1': 2.97500680936e-07,
  'acbyhs4.r2b2': -1.04121190725e-07,
  'acbyhs5.l2b1': -2.37517327511e-07,
  'acbyhs5.l2b2': -8.31782737488e-08},
 'on_a8': {'acbcv5.l8b2': 0.0,
  'acbcvs5.l8b1': -8.21760660637e-08,
  'acbcvs5.l8b2': -2.50223614299e-07,
  'acbyv4.l8b1': 0.0,
  'acbyv4.r8b2': 0.0,
  'acbyv5.r8b1': 0.0,
  'acbyvs4.l8b1': -3.28608937393e-07,
  'acbyvs4.l8b2': -4.88045166324e-09,
  'acbyvs4.r8b1': -2.49723205992e-08,
  'acbyvs4.r8b2': 3.29247846628e-07,
  'acbyvs5.r8b1': 2.35889032138e-07,
  'acbyvs5.r8b2': 8.24862621483e-08},
 'on_ov1': {'acbcv5.l1b1': 2.04334717035e-05,
  'acbcv5.r1b2': 1.77021721706e-05,
  'acbcv6.l1b2': 1.38877440621e-05,
  'acbcv6.r1b1': 1.38099768884e-05,
  'acbcv7.l1b1': -3.57270190667e-05,
  'acbcv7.r1b2': -3.43980837178e-05,
  'acbcv8.l1b2': -3.67951166259e-05,
  'acbcv8.r1b1': -3.54017584491e-05,
  'acbyvs4.l1b1': 2.18954414715e-05,
  'acbyvs4.l1b2': 1.48813814753e-05,
  'acbyvs4.r1b1': 1.47980502321e-05,
  'acbyvs4.r1b2': 1.89687234898e-05},
 'on_ov2': {'acbcv5.r2b1': -8.09924595773e-06,
  'acbcv6.l2b1': -6.04651281594e-06,
  'acbcv6.r2b2': -1.13045475282e-05,
  'acbcv7.l2b2': -1.99307633032e-05,
  'acbcv7.r2b1': -2.02127351779e-05,
  'acbcvs5.r2b1': -8.09924595773e-06,
  'acbcvs5.r2b2': -3.77186132907e-05,
  'acbyv4.l2b1': 6.47912745271e-06,
  'acbyv4.r2b2': 1.2113362935e-05,
  'acbyv5.l2b2': -2.35466128048e-06,
  'acbyvs4.l2b1': 6.47912745271e-06,
  'acbyvs4.l2b2': 4.16028721622e-05,
  'acbyvs4.r2b1': 4.21914517202e-05,
  'acbyvs4.r2b2': 1.2113362935e-05,
  'acbyvs5.l2b1': -4.64522709319e-05,
  'acbyvs5.l2b2': -2.35466128048e-06},
 'on_oh5': {'acbch5.l5b2': 2.34391469357e-05,
  'acbch5.r5b1': 2.09738881576e-05,
  'acbch6.l5b1': 1.00210256974e-05,
  'acbch6.r5b2': 1.12146824409e-05,
  'acbch7.l5b2': -3.61933012237e-05,
  'acbch7.r5b1': -3.89631993797e-05,
  'acbch8.l5b1': -2.19370052312e-05,
  'acbch8.r5b2': -2.35263402012e-05,
  'acbyhs4.l5b1': 1.07380079522e-05,
  'acbyhs4.l5b2': 2.51161661279e-05,
  'acbyhs4.r5b1': 2.24745235292e-05,
  'acbyhs4.r5b2': 1.20170681991e-05},
 'on_oh8': {'acbch5.l8b1': -3.75977649657e-06,
  'acbch6.l8b2': -1.18221036945e-05,
  'acbch6.r8b1': -7.34487460812e-06,
  'acbch7.l8b1': -2.14223976213e-05,
  'acbch7.r8b2': -1.76714732648e-05,
  'acbchs5.l8b1': -3.75977649657e-06,
  'acbchs5.l8b2': -2.59949391155e-05,
  'acbyh4.l8b2': 1.26679491018e-05,
  'acbyh4.r8b1': 7.87038416337e-06,
  'acbyh5.r8b2': -4.73283800199e-06,
  'acbyhs4.l8b1': 4.47164644971e-05,
  'acbyhs4.l8b2': 1.26679491018e-05,
  'acbyhs4.r8b1': 7.87038416337e-06,
  'acbyhs4.r8b2': 3.68868985081e-05,
  'acbyhs5.r8b1': -4.02874597858e-05,
  'acbyhs5.r8b2': -4.73283800199e-06}}


configs = {}

for nn in dct.keys():

    configs[nn] = {}

    purpose = ('sep' if '_sep' in nn
          else   'o' if '_o' in nn
          else 'a' if '_a' in nn
          else 'x' if '_x' in nn
          else None)
    assert purpose is not None

    configs[nn]['purpose'] = purpose

    # Determine the plane
    one_corrector = list(dct[nn].keys())[0]
    if 'h' in one_corrector and 'v' in one_corrector:
        raise ValueError(f'Cannot determine plane for {nn}')
    plane = 'x' if 'h' in one_corrector else 'y'
    configs[nn]['plane'] = plane

    # Determine ip
    tmp = one_corrector.split('.')[1][1]
    assert tmp in ['1', '2', '5', '8']
    ip = int(tmp)
    configs[nn]['ip'] = ip

