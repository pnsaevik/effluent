def rho(temp, salt, depth):
    """Simplified equation of state, from the NEMO documentation"""

    rho_0 = 1026         # seawater reference density
    a0 = 1.6550e-1       # linear thermal expansion coeff.
    b0 = 7.6554e-1       # linear haline expansion coeff.
    lambda1 = 5.9520e-2  # cabbeling coeff. in T^2
    lambda2 = 5.4914e-4  # cabbeling coeff. in S^2
    nu = 2.4341e-3       # cabbeling coeff. in T S
    mu1 = 1.4970e-4      # thermobaric coeff. in T
    mu2 = 1.1090e-5      # thermobaric coeff. in S

    t_a = temp - 10
    s_a = salt - 35
    z = depth

    da_times_rho0 = (
            -a0 * (1 + 0.5 * lambda1 * t_a + mu1 * z) * t_a
            + b0 * (1 - 0.5 * lambda2 * s_a - mu2 * z) * s_a
            - nu * t_a * s_a
    )

    return rho_0 + da_times_rho0
