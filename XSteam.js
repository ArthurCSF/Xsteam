/* **********************************************************************************************************
 ** Water and steam properties according to IAPWS IF-97                                                     *
 ** By Magnus Holmgren, www.x-eng.com                                                                       *
 ** The steam tables are free and provided as is.                                                           *
 ** We take no responsibilities for any errors in the code or damage thereby.                               *
 ** You are free to use, modify and distribute the code as long as authorship is properly acknowledged.     *
 ** Please notify me at magnus@x-eng.com if the code is used in commercial applications                     *
 ** The conversion of the code to JavaScript was performed by Arthur Colombo                                *
 ** Contact e-mail: arthur.colombo@outlook.com                                                              *
 ************************************************************************************************************
 **
 ** The code is also avalibale for matlab at www.x-eng.com
 **
 ** Contents.
 ** 1 Calling functions
 ** 1.1
 ** 1.2 Temperature (T)
 ** 1.3 Pressure (p)
 ** 1.4 Enthalpy (h)
 ** 1.5 Specific Volume (v)
 ** 1.6 Density (rho)
 ** 1.7 Specific entropy (s)
 ** 1.8 Specific internal energy (u)
 ** 1.9 Specific isobaric heat capacity (Cp)
 ** 1.10 Specific isochoric heat capacity (Cv)
 ** 1.11 Speed of sound
 ** 1.12 Viscosity
 ** 1.13 Prandtl
 ** 1.14 Kappa
 ** 1.15 Surface tension
 ** 1.16 Heat conductivity
 ** 1.17 Vapour fraction
 ** 1.18 Vapour Volume Fraction
 **
 ** 2 IAPWS IF 97 Calling functions
 ** 2.1 Functions for region 1
 ** 2.2 Functions for region 2
 ** 2.3 Functions for region 3
 ** 2.4 Functions for region 4
 ** 2.5 Functions for region 5
 **
 ** 3 Region Selection
 ** 3.1 Regions as a function of pT
 ** 3.2 Regions as a function of ph
 ** 3.3 Regions as a function of ps
 ** 3.4 Regions as a function of hs
 ** 3.5 Regions as a function of p and rho
 **
 ** 4 Region Borders
 ** 4.1 Boundary between region 1 and 3.
 ** 4.2 Region 3. pSat_h and pSat_s
 ** 4.3 Region boundary 1to3 and 3to2 as a functions of s
 **
 ** 5 Transport properties
 ** 5.1 Viscosity (IAPWS formulation 1985)
 ** 5.2 Thermal Conductivity (IAPWS formulation 1985)
 ** 5.3 Surface Tension
 **
 ** 6 Units
 */

/*************************************************************************************************************
 ** 1 Calling functions                                                                                      *
 ************************************************************************************************************/

//************************************************************************************************************
//** 1.1
//************************************************************************************************************

//************************************************************************************************************
//** 1.2 Temperature (T)

function Tsat_p(p) {
  p = toSIunit_p(p);

  if (p >= 0.000611657 && p <= 22.06395 + 0.001) {
    return fromSIunit_T(T4_p(p));
  } else {
    return Error("Invalid value");
  }
}

function Tsat_s(s) {
  s = toSIunit_s(s);

  if (s > -0.0001545495919 && s < 9.155759395) {
    return fromSIunit_T(T4_s(p4_s(s)));
  } else {
    return Error("Invalid value");
  }
}

function T_ph(p, h) {
  p = toSIunit_p(p);
  h = toSIunit_h(h);

  switch (region_ph(p, h)) {
    case 1:
      return fromSIunit_T(T1_ph(p, h));
    case 2:
      return fromSIunit_T(T2_ph(p, h));
    case 3:
      return fromSIunit_T(T3_ph(p, h));
    case 4:
      return fromSIunit_T(T4_p(p));
    case 5:
      return fromSIunit_T(T5_ph(p, h));
    default:
      return Error("Invalid value");
  }
}

function T_ps(p, s) {
  p = toSIunit_p(p);
  s = toSIunit_s(s);

  switch (region_ps(p, s)) {
    case 1:
      return fromSIunit_T(T1_ps(p, s));
    case 2:
      return fromSIunit_T(T2_ps(p, s));
    case 3:
      return fromSIunit_T(T3_ps(p, s));
    case 4:
      return fromSIunit_T(T4_p(p));
    case 5:
      return fromSIunit_T(T5_ps(p, s));
    default:
      return Error("Invalid value");
  }
}

function T_hs(h, s) {
  h = toSIunit_p(h);
  s = toSIunit_s(s);

  switch (region_hs(h, s)) {
    case 1:
      return fromSIunit_T(T1_ph(p1_hs(h, s), h));
    case 2:
      return fromSIunit_T(T2_ph(p2_hs(h, s), h));
    case 3:
      return fromSIunit_T(T3_ph(p3_hs(h, s), h));
    case 4:
      return fromSIunit_T(T4_p(p4_hs(h, s)));
    case 5:
      return Error("Invalid value"); // Functions of hs is not implemented in region 5
    default:
      return Error("Invalid value");
  }
}
//***********************************************************************************************************

//***********************************************************************************************************
//** 1.3 Pressure

function psat_T(T) {
  T = toSIunit_T(T);

  if (T <= 647.096 && T > 273.15) {
    return fromSIunit_p(p4_T(T));
  } else {
    return Error("Invalid value");
  }
}

function psat_s(s) {
  s = toSIunit_s(s);

  if (s > -0.0001545495919 && s < 9.155759395) {
    return fromSIunit_p(p4_s(s));
  } else {
    return Error("Invalid value");
  }
}

function p_hs(h, s) {
  h = toSIunit_h(h);
  s = toSIunit_s(s);

  switch (region_hs(h, s)) {
    case 1:
      return fromSIunit_p(p1_hs(h, s));
    case 2:
      return fromSIunit_p(p2_hs(h, s));
    case 3:
      return fromSIunit_p(p3_hs(h, s));
    case 4:
      return fromSIunit_p(p4_T(T4_hs(h, s)));
    case 5:
      return Error("Invalid value"); // Functions of hs is not implemented in region 5
    default:
      return Error("Invalid value");
  }
}

function p_hrho(h, rho) {
  // Not valid for water or supercritical since water rho does not change very much with p.
  // Uses iteration to find p.
  let High_Bound = fromSIunit_p(100);
  let Low_Bound = fromSIunit_p(0.000611657);
  let p = fromSIunit_p(10);
  let rhos = 1 / v_ph(p, h);

  while (Math.abs(rho - rhos) > 0.0000001) {
    rhos = 1 / v_ph(p, h);
    if (rhos >= rho) {
      High_Bound = p;
    } else {
      Low_Bound = p;
    }
    p = (Low_Bound + High_Bound) / 2;
  }
  return p;
}

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.4 Enthalpy (h)

function hV_p(p) {
  p = toSIunit_p(p);

  if (p > 0.000611657 && p < 22.06395) {
    return fromSIunit_h(h4V_p(p));
  } else {
    return Error("Invalid value");
  }
}

function hL_p(p) {
  p = toSIunit_p(p);

  if (p > 0.000611657 && p < 22.06395) {
    return fromSIunit_h(h4L_p(p));
  } else {
    return Error("Invalid value");
  }
}

function hV_T(T) {
  T = toSIunit_T(T);

  if (T > 273.15 && T < 647.096) {
    return fromSIunit_h(h4V_p(p4_T(T)));
  } else {
    return Error("Invalid value");
  }
}

function hL_T(T) {
  T = toSIunit_T(T);

  if (T > 273.15 && T < 647.096) {
    return fromSIunit_h(h4L_p(p4_T(T)));
  } else {
    return Error("Invalid value");
  }
}

function h_pT(p, T) {
  p = toSIunit_p(p);
  T = toSIunit_T(T);

  switch (region_pT(p, T)) {
    case 1:
      return fromSIunit_h(h1_pT(p, T));
    case 2:
      return fromSIunit_h(h2_pT(p, T));
    case 3:
      return fromSIunit_h(h3_pT(p, T));
    case 4:
      return Error("Invalid value");
    case 5:
      return fromSIunit_h(h5_pT(p, T));
    default:
      return Error("Invalid value");
  }
}

function h_ps(p, s) {
  p = toSIunit_p(p);
  s = toSIunit_s(s);

  switch (region_ps(p, s)) {
    case 1:
      return fromSIunit_h(h1_pT(p, T1_ps(p, s)));
    case 2:
      return fromSIunit_h(h2_pT(p, T2_ps(p, s)));
    case 3:
      return fromSIunit_h(h3_rhoT(1 / v3_ps(p, s), T3_ps(p, s)));
    case 4:
      const xs = x4_ps(p, s);
      return fromSIunit_h(xs * h4V_p(p) + (1 - xs) * h4L_p(p));
    case 5:
      return fromSIunit_h(h5_pT(p, T5_ps(p, s)));
    default:
      return Error("Invalid value");
  }
}

function h_px(p, x) {
  p = toSIunit_p(p);
  x = toSIunit_x(x);

  if ((x > 1) | (x < 0) | (p >= 22.064)) {
    return Error("Invalid value");
  }
  const hL = h4L_p(p);
  const hV = h4V_p(p);
  return fromSIunit_h(hL + x * (hV - hL));
}

function h_Tx(T, x) {
  T = toSIunit_T(T);
  x = toSIunit_x(x);

  if ((x > 1) | (x < 0) | (T >= 647.096)) {
    return Error("Invalid value");
  }

  const p = p4_T(T);
  const hL = h4L_p(p);
  const hV = h4V_p(p);
  return fromSIunit_h(hl + x * (hV - hL));
}

function h_prho(p, rho) {
  p = toSIunit_p(p);
  rho = 1 / toSIunit_v(1 / rho);

  switch (region_prho(p, rho)) {
    case 1:
      return fromSIunit_h(h1_pT(p, T1_prho(p, rho)));
    case 2:
      return fromSIunit_h(h2_pT(p, T2_prho(p, rho)));
    case 3:
      return fromSIunit_h(h3_rhoT(rho, T3_prho(p, rho)));
    case 4:
      let vV = 0;
      let vL = 0;
      const hV = h4V_p(p);
      const hL = h4L_p(p);
      if (p < 16.529) {
        vV = v2_pT(p, T4_p(p));
        vL = v1_pT(p, T4_p(p));
      } else {
        vV = v3_ph(p, h4V_p(p));
        vL = v3_ph(p, h4L_p(p));
      }
      const x = (1 / rho - vL) / (vV - vL);
      return fromSIunit_h((1 - x) * hL + x * hV);
    case 5:
      return fromSIunit_h(h5_pT(p, T5_prho(p, rho)));
    default:
      return Error("Invalid value");
  }
}

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.5 Specific Volume (v)

function vV_p(p) {
  p = toSIunit_p(p);

  if (p > 0.00611657 && p < 22.06395) {
    if (p < 16.529) {
      return fromSIunit_v(v2_pT(p, T4_p(p)));
    } else {
      return fromSIunit_v(v3_ph(p, h4V_p(p)));
    }
  } else {
    return Error("Invalid value");
  }
}

function vL_p(p) {
  p = toSIunit_p(p);

  if (p > 0.00611657 && p < 22.06395) {
    if (p < 16.529) {
      return fromSIunit_v(v2_pT(p, T4_p(p)));
    } else {
      return fromSIunit_v(v3_ph(p, h4L_p(p)));
    }
  } else {
    return Error("Invalid value");
  }
}

function vV_T(T) {
  T = toSIunit_T(T);

  if (T > 273.15 && T < 647.096) {
    if (T <= 623.15) {
      return fromSIunit_v(v2_pT(p4_T(T), T));
    } else {
      return fromSIunit_v(v3_ph(p4_T(T), h4V_p(p4_T(T))));
    }
  } else {
    return Error("Invalid value");
  }
}

function vL_T(T) {
  T = toSIunit_T(T);

  if (T > 273.15 && T < 647.096) {
    if (T <= 623.15) {
      return fromSIunit_v(v1_pT(p4_T(T), T));
    } else {
      return fromSIunit_v(v3_ph(p4_T(T), h4L_p(p4_T(T))));
    }
  } else {
    return Error("Invalid value");
  }
}

function v_pT(p, T) {
  p = toSIunit_p(p);
  T = toSIunit_T(T);

  switch (region_pT(p, T)) {
    case 1:
      return fromSIunit_v(v1_pT(p, T));
    case 2:
      return fromSIunit_v(v2_pT(p, T));
    case 3:
      return fromSIunit_v(v3_ph(p, h3_pT(p, T)));
    case 4:
      return Error("Invalid value");
    case 5:
      fromSIunit_v(v5_pT(p, T));
      return;
    default:
      return Error("Invalid value");
  }
}

function v_ph(p, h) {
  p = toSIunit_p(p);
  h = toSIunit_h(h);

  switch (region_ph(p, h)) {
    case 1:
      return fromSIunit_v(v1_pT(p, T1_ph(p, h)));
    case 2:
      return fromSIunit_v(v2_pT(p, T2_ph(p, h)));
    case 3:
      return fromSIunit_v(v3_ph(p, h));
    case 4:
      const xs = x4_ph(p, h);
      let v4V;
      let v4L;
      if (p < 16.529) {
        v4V = v2_pT(p, T4_p(p));
        v4L = v1_pT(p, T4_p(p));
      } else {
        v4V = v3_ph(p, h4V_p(p));
        v4L = v3_ph(p, h4L_p(p));
      }
      return fromSIunit_v(xs * v4V + (1 - xs) * v4L);
    case 5:
      return fromSIunit_v(v5_pT(p, T5_ph(p, h)));
    default:
      return Error("Invalid value");
  }
}

function v_ps(p, s) {
  p = toSIunit_p(p);
  s = toSIunit_s(s);

  switch (region_ps(p, s)) {
    case 1:
      return fromSIunit_v(v1_pT(p, T1_ps(p, s)));
    case 2:
      return fromSIunit_v(v2_pT(p, T2_ps(p, s)));
    case 3:
      return fromSIunit_v(v3_ps(p, s));
    case 4:
      const xs = x4_ps(p, s);
      let v4V;
      let v4L;
      if (p < 16.529) {
        v4V = v2_pT(p, T4_p(p));
        v4L = v1_pT(p, T4_p(p));
      } else {
        v4V = v3_ph(p, h4V_p(p));
        v4L = v3_ph(p, h4L_p(p));
      }
      return fromSIunit_v(xs * v4V + (1 - xs) * v4L);
    case 5:
      return fromSIunit_v(v5_pT(p, T5_ps(p, s)));
    default:
      return Error("Invalid value");
  }
}

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.6 Density (rho)
//  Density is calculated as 1/v

function rhoV_p(p) {
  return 1 / vV_p(p);
}

function rhoL_p(p) {
  return 1 / vL_p(p);
}

function rhoL_T(T) {
  return 1 / vL_T(T);
}

function rhoV_T(T) {
  return 1 / vV_T(T);
}

function rho_pT(p, T) {
  return 1 / v_pT(p, T);
}

function rho_ph(p, h) {
  return 1 / v_ph(p, h);
}

function rho_ps(p, s) {
  return 1 / v_ps(p, s);
}

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.7 Specific entropy (s)

function sV_p(p) {
  p = toSIunit_p(p);

  if (p > 0.00611657 && p < 22.06395) {
    if (p < 16.529) {
      return fromSIunit_s(s2_pT(p, T4_p(p)));
    } else {
      return fromSIunit_s(s3_rhoT(1 / v3_ph(p, h4V_p(p)), T4_p(p)));
    }
  } else {
    return Error("Invalid value");
  }
}

function sL_p(p) {
  p = toSIunit_p(p);

  if (p > 0.00611657 && p < 22.06395) {
    if (p < 16.529) {
      return fromSIunit_s(s1_pT(p, T4_p(p)));
    } else {
      return fromSIunit_s(s3_rhoT(1 / v3_ph(p, h4L_p(p)), T4_p(p)));
    }
  } else {
    return Error("Invalid value");
  }
}

function sV_T(T) {
  T = toSIunit_T(T);

  if (T > 273.15 && T < 647.096) {
    if (T <= 623.15) {
      return fromSIunit_s(s2_pT(p4_T(T), T));
    } else {
      fromSIunit_s(s3_rhoT(1 / v3_ph(p4_T(T), h4V_p(p4_T(T))), T));
    }
  } else {
    return Error("Invalid value");
  }
}

function sL_T(T) {
  T = toSIunit_T(T);

  if (T > 273.15 && T < 647.096) {
    if (T <= 623.15) {
      return fromSIunit_s(s1_pT(p4_T(T), T));
    } else {
      fromSIunit_s(s3_rhoT(1 / v3_ph(p4_T(T), h4L_p(p4_T(T))), T));
    }
  } else {
    return Error("Invalid value");
  }
}

function s_pT(p, T) {
  p = toSIunit_p(p);
  T = toSIunit_T(T);

  switch (region_pT(p, T)) {
    case 1:
      return fromSIunit_s(s1_pT(p, T));
    case 2:
      return fromSIunit_s(s2_pT(p, T));
    case 3:
      return fromSIunit_s(s3_rhoT(1 / v3_ph(p, h3_pT(p, T)), T));
    case 4:
      return Error("Invalid value");
    case 5:
      return fromSIunit_s(s5_pT(p, T));
    default:
      return Error("Invalid value");
  }
}

function s_ph(p, h) {
  p = toSIunit_p(p);
  h = toSIunit_h(h);

  switch (region_ph(p, h)) {
    case 1:
      return fromSIunit_s(s1_pT(p, T1_ph(p, h)));
    case 2:
      return fromSIunit_s(s2_pT(p, T2_ph(p, h)));
    case 3:
      return fromSIunit_s(s3_rhoT(1 / v3_ph(p, h), T3_ph(p, h)));
    case 4:
      const Ts = T4_p(p);
      const xs = x4_ph(p, h);
      let v4V;
      let s4V;
      let v4L;
      let s4L;
      if(p < 16.529){
        s4V = s2_pT(p, Ts);
        s4L = s1_pT(p, Ts);
      }else{
        v4V = v3_ph(p, h4V_p(p));
        s4V = s3_rhoT(1 / v4V, Ts);
        v4L = v3_ph(p, h4L_p(p));
        s4L = s3_rhoT(1 / v4L, Ts);
      }
      return fromSIunit_s((xs * s4V + (1 - xs) * s4L))
    case 5:
      return fromSIunit_s(s5_pT(p, T5_ph(p, h)));
    default:
      return Error("Invalid value");
  }
}
//***********************************************************************************************************

//***********************************************************************************************************
//** 1.8 Specific internal energy (u)



//***********************************************************************************************************

//***********************************************************************************************************
//** 1.9 Specific isobaric heat capacity (Cp)

function CpV_p(p){
  p= toSIunit_p(p);
  if(p > 0.000611657 && p < 22.06395){
    if(p < 16.529){
      return fromSIunit_Cp(Cp2_pT(p, T4_p(p)));
    } else {
      return fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p, h4V_p(p))),T4_p(p)));
    } 
  } else {
    return Error("Invalid Value");
  }
}

function CpL_p(p){
  p = toSIunit_p(p);
  if(p > 0.000611657 && p < 22.06395){
    if(p < 16.529){
      return fromSIunit_Cp(Cp1_pT(p, T4_p(p)));
    } else {
      const T = T4_p(p);
      const h = h4L_p(p);
      const v = v3_ph(p, h4L_p(p));
      return fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p)));
    }
  } else {
    return Error("Invalid Value");
  }
}

function CpV_T(T){
T = toSIunit_T(T);
if(T > 273.15 && T < 647.096){
  if(T <= 623.15){
    return fromSIunit_Cp(Cp2_pT(p4_T(T), T));
} else {
  return fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p4_T(T), h4V_p(p4_T(T)))), T));
}
} else {
  return Error("Invalid Value");
}
}


//***********************************************************************************************************

//***********************************************************************************************************
//** 1.10 Specific isochoric heat capacity (Cv)

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.11 Speed of sound

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.12 Viscosity

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.13 Prandtl

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.14 Kappa

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.15 Surface tension

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.16 Heat conductivity

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.17 Vapour fraction

//***********************************************************************************************************

//***********************************************************************************************************
//** 1.18 Vapour Volume Fraction

//***********************************************************************************************************

/************************************************************************************************************
 ** 2 IAPWS IF 97 Calling functions                                                                          *
 ************************************************************************************************************/

//***********************************************************************************************************
//** 2.1 Functions for region 1

//***********************************************************************************************************

//***********************************************************************************************************
//** 2.2 Functions for region 2

//***********************************************************************************************************

//***********************************************************************************************************
//** 2.3 Functions for region 3

//***********************************************************************************************************

//***********************************************************************************************************
//** 2.4 Functions for region 4

//***********************************************************************************************************

//***********************************************************************************************************
//** 2.5 Functions for region 5

//***********************************************************************************************************

/************************************************************************************************************
 ** 3 Region Selection                                                                                       *
 ************************************************************************************************************/

//***********************************************************************************************************
//** 3.1 Regions as a function of pT

//***********************************************************************************************************

//***********************************************************************************************************
//** 3.2 Regions as a function of ph

//***********************************************************************************************************

//***********************************************************************************************************
//** 3.3 Regions as a function of ps

//***********************************************************************************************************

//***********************************************************************************************************
//** 3.4 Regions as a function of hs

//***********************************************************************************************************

//***********************************************************************************************************
//** 3.5 Regions as a function of p and rho

//***********************************************************************************************************

/************************************************************************************************************
 ** 4 Region Borders                                                                                         *
 ************************************************************************************************************/

//***********************************************************************************************************
//** 4.1 Boundary between region 1 and 3

//***********************************************************************************************************

//***********************************************************************************************************
//** 4.2 Region 3. pSat_h and pSat_s

//***********************************************************************************************************

//***********************************************************************************************************
//** 4.3 Region boundary 1to3 and 3to2 as a functions of s

//***********************************************************************************************************

/************************************************************************************************************
 ** 5 Transport properties                                                                                   *
 ************************************************************************************************************/

//***********************************************************************************************************
//** 5.1 Viscosity (IAPWS formulation 1985)

//***********************************************************************************************************

//***********************************************************************************************************
//** 5.2 Thermal Conductivity (IAPWS formulation 1985)

//***********************************************************************************************************

//***********************************************************************************************************
//** 5.3 Surface Tension

//***********************************************************************************************************

/************************************************************************************************************
 ** 6 Units                                                                                                  *
 ************************************************************************************************************/
