module PythonToolsExt

using PythonCall, SunnyAnalysisTools


struct TAVISpec <: SunnyAnalysisTools.AbstractInstrumentSpec
    name             :: String
    tavi_instrument  :: Py
    params           :: Dict{String, Any}
end


function SunnyAnalysisTools.TAVISpec(name, instrumentfile::String; Ef, samplefile = nothing, params=Dict())
    TAS = pyimport("tavi.instrument.tas").TAS
    tavi_instrument = TAS(fixed_ef=Ef)
    tavi_instrument.load_instrument_params_from_json(instrumentfile)
    if !isnothing(samplefile)
        Sample = pyimport("tavi.sample").Sample 
        sample = Sample.from_json(samplefile)
        tavi_instrument.mount_sample(sample)
    end
    return TAVISpec(name, tavi_instrument, params)
end

function SunnyAnalysisTools.TripleAxisMC(path::TripleAxisPath, instrument::TAVISpec; N)
    (; HKLs, Es, projection) = path
    (; tavi_instrument) = instrument

    projection_py = tuple([tuple(col...) for col in eachcol(projection)]...)
    resmats_py = [tavi_instrument.cooper_nathans(hkl=pylist([tuple(HKL...)]), en=pylist([E]), projection=projection_py) for (HKL, E) in zip(HKLs, Es)]
    Ks = [pyconvert(Array{Float64, 2}, res_elipse.mat) for res_elipse in resmats_py]
    println(length(Ks))

    SunnyAnalysisTools.TripleAxisMC{1}(path, Ks, N)
end

# Assuming a 2D contour plot is desired (take product of Q and E points).
function SunnyAnalysisTools.TripleAxisMC(path::TripleAxis2DContour, instrument::TAVISpec; N)
    (; HKLs, Es, projection) = path
    (; tavi_instrument) = instrument

    HKLs_py = pylist([tuple(hkl...) for hkl in HKLs])
    Es_py = pylist(Es) 
    projection_py = tuple([tuple(col...) for col in eachcol(projection)]...)

    rez_list = tavi_instrument.cooper_nathans(hkl=HKLs_py, en=Es_py, projection=projection_py)

    Ks = [pyconvert(Array{Float64, 2}, res_elipse.mat) for res_elipse in rez_list]
    hkl = [pyconvert(Vector, res_elipse.hkl) for res_elipse in rez_list]
    Es = [pyconvert(Float64, res_elipse.en) for res_elipse in rez_list]

    SunnyAnalysisTools.TripleAxisMC{2}(path, Ks, N)
end

function SunnyAnalysisTools.TripleAxisGrid(path::TripleAxis2DContour, instrument::TAVISpec; nsigmas=3, counts=3*ones(3))
    (; HKLs, Es, projection) = path
    (; tavi_instrument) = instrument

    HKLs_py = pylist([tuple(hkl...) for hkl in HKLs])
    Es_py = pylist(Es) 
    projection_py = tuple([tuple(col...) for col in eachcol(projection)]...)

    rez_list = tavi_instrument.cooper_nathans(hkl=HKLs_py, en=Es_py, projection=projection_py)

    Ks = [pyconvert(Array{Float64, 2}, res_elipse.mat) for res_elipse in rez_list]
    hkl = [pyconvert(Vector, res_elipse.hkl) for res_elipse in rez_list]
    Es = [pyconvert(Float64, res_elipse.en) for res_elipse in rez_list]

    if length(counts) == 1
        counts = ones(3)*counts
    end
    SunnyAnalysisTools.TripleAxisGrid{2}(path, Ks, nsigmas, Tuple(counts))
end

function SunnyAnalysisTools.cncs(; Ei, Δθ = 1.5)
    Instruments = pyimport("PyChop.Instruments")
    CNCS = Instruments.Instrument("CNCS")

    Δθ = Δθ * (π/180) # Convert to radians
    distances = CNCS.chopper_system.getDistances()
    x0, _, x1, x2, xm = map(val -> pyconvert(Float64, val), distances)  # xa not used
    tsqmod = CNCS.moderator.getWidthSquared(Ei) |> x -> pyconvert(Float64, x)
    chopper_width = CNCS.chopper_system.getWidthSquared(Ei)[0] |> x -> pyconvert(Float64, x) |> sqrt

    L1, L2, L3 = x0 - xm, x1, x2

    # Rescale the moderator width to its effective value at the first (pulse-shaping)
    # chopper, then take the tighter of that and the pulse-shaping chopper's own
    # opening-time width -- mirrors PyChop's own Instrument.getVanVar (Instruments.py)
    # rather than always assuming the rescaled moderator width is binding.
    frac_dist = 1 - (xm/x0)
    tsmeff = tsqmod * frac_dist^2
    tsqchp1 = pyconvert(Union{Float64,Nothing}, CNCS.chopper_system.getWidthSquared(Ei)[1])
    tsqp = isnothing(tsqchp1) ? tsmeff : (tsqchp1 > tsmeff ? tsmeff : tsqchp1)
    Δtp = sqrt(tsqp)
    Δtc = chopper_width
    # Detector.getWidth returns a length (m, the He-3 tube absorption-depth spread), not
    # directly a time -- PyChop's own getVanVar converts via (1/vf)^2, so divide by the
    # elastic-line final velocity to get an actual time width comparable to Δtp/Δtc.
    detector_width_m = CNCS.detector.getWidth(Ei) |> x -> pyconvert(Float64, x)
    vf = SunnyAnalysisTools.energy_to_velocity(Ei)
    Δtd = detector_width_m / vf

    return SunnyAnalysisTools.ChopperSpec("CNCS", Ei, L1, L2, L3, Δtp, Δtc, Δtd, Δθ)
end

function SunnyAnalysisTools.spins()
    params = Dict{String, Any}()
    return SunnyAnalysisTools.TripleAxisSpec("SPINS", params)
end

function SunnyAnalysisTools.hmi(; hires=false)
    params = Dict{String, Any}()
    if hires
        params["resolution"] = "high"
    else
        params["resolution"] = "low"
    end
    return SunnyAnalysisTools.TripleAxisSpec("HMI", params)
end



end