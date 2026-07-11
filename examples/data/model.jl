function bamno_model(crystal;
    J0 = 1.642,
    J1 = 0.118,
    J2 = 0.256,
    J3 = 0.142,
    J4 = 0.037,
    D  = -0.032,
    dims = (1,1,1),
)
    # Define the system and interactions
    sys_origin = System(crystal, [1 => Moment(; s=1, g=2)], :SUN; dims)
    set_exchange!(sys_origin, J0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys_origin, J1, Bond(2, 3, [0, 0, 0]))
    set_exchange!(sys_origin, J2, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys_origin, J3, Bond(1, 2, [1, 0, 0]))
    set_exchange!(sys_origin, J4, Bond(4, 5, [0, 1, 0]))
    set_onsite_coupling!(sys_origin, S -> D*S[3]^2, 1)

    # Entangle relevant sites
    dimers = [(1, 2), (3, 4), (5, 6)]
    sys = Sunny.EntangledSystem(sys_origin, dimers)

    # Initialize in singlet state
    singlet = [0.0 + 0.0im, 0.0 + 0.0im, -0.5773502691896271 + 0.0im, 0.0 + 0.0im, 0.5773502691896258 + 0.0im, 0.0 + 0.0im, -0.5773502691896244 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]
    for unit in Sunny.eachunit(sys)
        set_coherent!(sys, singlet, unit)
    end

    return (; sys, crystal)
end

# Make a SunnyAnalysisTools model
function make_swt_model(crystal, params)
    (; sys) = bamno_model(crystal; params...)
    formfactors = [1 => FormFactor("Mn2")]
    measure = ssf_perp(sys; formfactors)
    minimize_energy!(sys)
    swt = SpinWaveTheory(sys; measure)
    return swt
end