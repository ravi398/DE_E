      double f(double x) {
      // Omega_ede(a) taken from eq. (10) in 1706.00730
       Omega_ede = (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(x,-3.*pba->w0_fld)))/(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(x,3.*pba->w0_fld))+ pba->Omega_EDE*(1.-pow(x,-3.*pba->w0_fld));

       // d Omega_ede / d a taken analytically from the above
       dOmega_ede_over_da = - pba->Omega_EDE* 3.*pba->w0_fld*pow(x,-3.*pba->w0_fld-1.)/(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(x,3.*pba->w0_fld))- (pba->Omega0_fld - pba->Omega_EDE*(1.-            pow(x,-3.*pba->w0_fld)))*(1.-pba->Omega0_fld)*3.*pba->w0_fld*pow(x,3.*pba->w0_fld-1.)/pow(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(x,3.*pba->w0_fld),2)
      + pba->Omega_EDE*3.*pba->w0_fld*pow(x,-3.*pba->w0_fld-1.);

         // find a_equality (needed because EDE tracks first radiation, then matter)
         Omega_r = pba->Omega0_g * (1. + 3.044 * 7./8.*pow(4./11.,4./3.)); // assumes LambdaCDM + eventually massive neutrinos so light that they are relativistic at equality; needs to be generalised later on.
        Omega_m = pba->Omega0_b;
        if (pba->has_cdm == _TRUE_) Omega_m += pba->Omega0_cdm;
        if (pba->has_idm == _TRUE_) Omega_m += pba->Omega0_idm;
        if (pba->has_dcdm == _TRUE_)
         class_stop(pba->error_message,"Early Dark Energy not compatible with decaying Dark Matter because we omitted to code the calculation of a_eq in that case, but it would not be difficult to add it if necessary, should be a matter of 5 minutes");
        a_eq = Omega_r/Omega_m; // assumes a flat universe with a=1 today

        // w_ede(a) taken from eq. (11) in 1706.00730
        return - dOmega_ede_over_da*x/Omega_ede/3./(1.-Omega_ede)+a_eq/3./(x+a_eq);
       }
