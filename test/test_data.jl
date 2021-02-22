sv_constants = GPSL1Constants()

sv1_data = GPSData(
    integrity_status_flag= false,
    TOW = 34810,
    alert_flag = false,
    anti_spoof_flag = true,
    trans_week = 67,
    codeonl2 = 1,
    ura = 2.0,
    svhealth = "000000",
    IODC = "0000011011",
    l2pcode = false,
    T_GD = -1.1175870895385742e-8,
    t_oc = 215984,
    a_f2 = 0.0,
    a_f1 = -8.640199666842818e-12,
    a_f0 = -0.0003322265110909939,
    IODE_Sub_2 = "00011011",
    C_rs = -15.6875,
    Δn = 5.142714214689423e-9,
    M_0 = -1.2713497555403992,
    C_uc = -6.761401891708374e-7,
    e = 0.014137965277768672,
    C_us= 4.753470420837402e-6,
    sqrt_A= 5153.64093208313,
    t_oe = 215984,
    fit_interval = false,
    AODO =13,
    C_ic =  3.129243850708008e-7,
    Ω_0 = -1.0654053987541545,
    C_is = 4.0978193283081055e-8,
    i_0 = 0.9525395470841136,
    C_rc = 276.71875,
    ω = -2.378063276589369,
    Ω_dot = -8.233200088688144e-9,
    IODE_Sub_3 = "00011011",
    IDOT = -1.2286226056345314e-10
)
sv1_decoder_struct = GNSSDecoderState(
    PRN = 7,
    buffer = BitArray(undef, 1502),
    data = sv1_data,
    constants = sv_constants,
    preamble_found =false,
    subframe_count =0,
    prev_30=false,
    prev_29=false,
    data_integrity=true,
    new_data_needed= false,
    subframes_decoded = [1,1,1,1,1],
    nb_prev = 0,
    num_bits_buffered = 187
)
sv1_code_phase = 18126.41231251371
sv1_carrier_phase = -0.019570794887840748
sv1_struct = SatelliteState(decoder_state = sv1_decoder_struct, code_phase = sv1_code_phase, carrier_phase = sv1_carrier_phase)



sv2_data = GPSData(
    integrity_status_flag = false,
    TOW =  34810,
    alert_flag =  false,
    anti_spoof_flag =  true,
    trans_week =  67,
    codeonl2 =  1,
    ura =  2.0,
    svhealth =  "000000",
    IODC =  "0001011001",
    l2pcode =  false,
    T_GD =  5.122274160385132e-9,
    t_oc =  216000,
    a_f2 =  0.0,
    a_f1 =  -1.3642420526593924e-12,
    a_f0 =  -4.171021282672882e-5,
    IODE_Sub_2 =  "01011001",
    C_rs =  148.90625,
    Δn =  4.43518474324698e-9,
    M_0 =  1.74993761786873,
    C_uc =  7.92182981967926e-6,
    e =  0.0053556435741484165,
    C_us =  5.660578608512878e-6,
    sqrt_A =  5153.701900482178,
    t_oe =  216000,
    fit_interval =  false,
    AODO =  31,
    C_ic =  1.0617077350616455e-7,
    Ω_0 =  1.0019742934168678,
    C_is =  8.940696716308594e-8,
    i_0 =  0.9693796725429337,
    C_rc =  272.0,
    ω =  -0.060505953824413546,
    Ω_dot =  -8.216770832915124e-9,
    IODE_Sub_3 =  "01011001",
    IDOT =  2.1858053332800384e-10
)
sv2_decoder_struct = GNSSDecoderState(
    PRN= 8,
    buffer = BitArray(undef, 1502),
    data = sv2_data,
    constants = sv_constants,
    preamble_found = false,
    subframe_count =  0,
    prev_30 =  false,
    prev_29 =  false,
    data_integrity =  true,
    new_data_needed =  false,
    subframes_decoded = [1,1,1,1,1],
    nb_prev =  1,
    num_bits_buffered =  188
)
sv2_code_phase = 1731.7849623335098
sv2_carrier_phase = 0.3726482419297099
sv2_struct = SatelliteState(decoder_state = sv2_decoder_struct, code_phase = sv2_code_phase, carrier_phase = sv2_carrier_phase)


sv3_data = GPSData(
    integrity_status_flag =  false,
    TOW =  34810,
    alert_flag =  false,
    anti_spoof_flag =  true,
    trans_week =  67,
    codeonl2 =  1,
    ura =  2.0,
    svhealth =  "000000",
    IODC =  "0000001001",
    l2pcode =  false,
    T_GD =  2.3283064365386963e-9,
    t_oc =  215984,
    a_f2 =  0.0,
    a_f1 =  -1.057287590811029e-11,
    a_f0 =  -0.00040573813021183014,
    IODE_Sub_2 =  "00001001",
    C_rs =  -128.0,
    Δn =  3.82944522605042e-9,
    M_0 =  -2.5953206787807925,
    C_uc =  -6.614252924919128e-6,
    e =  0.005705653806217015,
    C_us =  1.1416152119636536e-5,
    sqrt_A =  5153.6776695251465,
    t_oe =  215984,
    fit_interval =  false,
    AODO =  11,
    C_ic =  -4.842877388000488e-8,
    Ω_0 =  3.1100423470271776,
    C_is =  -3.166496753692627e-8,
    i_0 =  0.9666236109533345,
    C_rc =  160.65625,
    ω =  -2.6449072365113864,
    Ω_dot =  -7.597102164084918e-9,
    IODE_Sub_3 =  "00001001",
    IDOT =  4.071598169835366e-11
)
sv3_decoder_struct = GNSSDecoderState(
    PRN =  10,
    buffer = BitArray(undef, 1502),
    data = sv3_data,
    constants = sv_constants,
    preamble_found =  false,
    subframe_count =  0,
    prev_30 =  true,
    prev_29 =  true,
    data_integrity =  true,
    new_data_needed =  false,
    subframes_decoded = [1,1,1,1,1],
    nb_prev =  1,
    num_bits_buffered =  188
)
sv3_code_phase = 2975.5753833700155
sv3_carrier_phase = 0.007562681101262569
sv3_struct = SatelliteState(decoder_state = sv3_decoder_struct, code_phase = sv3_code_phase, carrier_phase = sv3_carrier_phase)


sv4_data = GPSData(
    integrity_status_flag =  false,
  TOW =  34810,
  alert_flag =  false,
  anti_spoof_flag =  true,
  trans_week =  67,
  codeonl2 =  1,
  ura =  2.0,
  svhealth =  "000000",
  IODC =  "0000010100",
  l2pcode =  false,
  T_GD =  -1.0710209608078003e-8,
  t_oc =  216000,
  a_f2 =  0.0,
  a_f1 =  2.6147972675971687e-12,
  a_f0 =  -0.00021605053916573524,
  IODE_Sub_2 =  "00010100",
  C_rs =  26.09375,
  Δn =  5.252004481353425e-9,
  M_0 =  0.768362779336754,
  C_uc =  1.4789402484893799e-6,
  e =  0.012468118453398347,
  C_us =  9.216368198394775e-6,
  sqrt_A =  5153.7172775268555,
  t_oe =  216000,
  fit_interval =  false,
  AODO =  28,
  C_ic =  -2.3469328880310059e-7,
  Ω_0 =  -2.2306191090736784,
  C_is =  6.891787052154541e-8,
  i_0 =  0.9282541736265442,
  C_rc =  186.90625,
  ω =  0.8934476147650858,
  Ω_dot =  -8.434994208508933e-9,
  IODE_Sub_3 =  "00010100",
  IDOT =  5.610948004220491e-10
)
sv4_decoder_struct = GNSSDecoderState(
    PRN =  15,
    buffer = BitArray(undef, 1502),
    data = sv4_data,
    constants = sv_constants,
    preamble_found =  false,
    subframe_count =  0,
    prev_30 =  false,
    prev_29 =  false,
    data_integrity =  true,
    new_data_needed =  false,
    subframes_decoded = [1,1,1,1,1],
    nb_prev =  0,
    num_bits_buffered =  187
)
sv4_code_phase = 17790.853982632394
sv4_carrier_phase = 0.2946765711531043
sv4_struct = SatelliteState(decoder_state = sv4_decoder_struct, code_phase = sv4_code_phase, carrier_phase = sv4_carrier_phase)


sv5_data = GPSData(
    integrity_status_flag = false,
    TOW = 34810,
    alert_flag = false,
    anti_spoof_flag = true,
    trans_week = 67,
    codeonl2 = 1,
    ura = 2.0,
    svhealth = "000000",
    IODC = "1011101101",
    l2pcode = false,
    T_GD = -7.916241884231567e-9,
    t_oc = 216000,
    a_f2 = 0.0,
    a_f1 = 9.208633855450898e-12,
    a_f0 = 0.0002515711821615696,
    IODE_Sub_2 = "11101101",
    C_rs = -81.9375,
    Δn = 4.641621913612317e-9,
    M_0 = 0.01827858501071997,
    C_uc = -4.407018423080444e-6,
    e = 0.0007737367413938046,
    C_us = 1.044943928718567e-6,
    sqrt_A = 5153.718004226685,
    t_oe = 216000,
    fit_interval = false,
    AODO = 31,
    C_ic = -2.60770320892334e-8,
    Ω_0 = 2.089468318582403,
    C_is = -2.9802322387695312e-8,
    i_0 = 0.9643063355495126,
    C_rc = 363.78125,
    ω = 2.6017948651652136,
    Ω_dot = -8.454637883889717e-9,
    IODE_Sub_3 = "11101101",
    IDOT = -3.8001582918463414e-10
)
sv5_decoder_struct = GNSSDecoderState(
    PRN =  18,
    buffer = BitArray(undef, 1502),
    data = sv5_data,
    constants = sv_constants,
    preamble_found= false,
    subframe_count= 0,
    prev_30 = false,
    prev_29 = true,
    data_integrity = true,
    new_data_needed = false,
    subframes_decoded = [1,1,1,1,1],
    nb_prev = 0,
    num_bits_buffered = 188
)

sv5_code_phase =  10736.222919350244
sv5_carrier_phase = 0.23463375493884087

sv5_struct = SatelliteState(decoder_state = sv5_decoder_struct, code_phase = sv5_code_phase, carrier_phase = sv5_carrier_phase)

satellite_states = [sv1_struct, sv2_struct, sv3_struct, sv4_struct, sv5_struct]
#test_dcs = Vector([sv1_struct, sv2_struct, sv3_struct, sv4_struct, sv5_struct])
#test_cops = Vector([sv1_code_phase, sv2_code_phase, sv3_code_phase, sv4_code_phase, sv5_code_phase])
#test_caps = Vector([sv1_carrier_phase, sv2_carrier_phase, sv3_carrier_phase, sv4_carrier_phase, sv5_carrier_phase])


