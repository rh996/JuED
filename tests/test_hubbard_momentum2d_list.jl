include("../src/EDMain.jl")
using PythonCall
plt = pyimport("matplotlib.pyplot")
using HDF5

using .EDMod

function test(nalpha, nbeta, Nkx, Nky)

    U = 5.0
    t = -1.0
    Eri = zeros(Float64, Nkx * Nky * 2, Nkx * Nky * 2, Nkx * Nky * 2, Nkx * Nky * 2)
    Onebody = zeros(Float64, Nkx * Nky * 2, Nkx * Nky * 2)

    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            for ikpx in 0:Nkx-1
                for ikpy in 0:Nky-1
                    for iqx in 0:Nkx-1
                        for iqy in 0:Nky-1
                            s1 = (mod(ikx + iqx, Nkx) + mod(iky + iqy, Nky) * Nkx) * 2
                            s2 = (ikx + iky * Nkx) * 2
                            s3 = (mod(ikpx - iqx, Nkx) + mod(ikpy - iqy, Nky) * Nkx) * 2 + 1
                            s4 = (ikpx + ikpy * Nkx) * 2 + 1
                            Eri[s1+1, s3+1, s2+1, s4+1] += U / (Nkx * Nky)
                        end
                    end
                end
            end
        end

    end


    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            ik = ikx + iky * Nkx
            Onebody[2*ik+1, 2*ik+1] = 2 * t * (cos(2 * pi * ikx / Nkx) + cos(2 * pi * iky / Nky))
            Onebody[2*ik+2, 2*ik+2] = 2 * t * (cos(2 * pi * ikx / Nkx) + cos(2 * pi * iky / Nky))
        end
    end


    modelparams = InputModel(nalpha, nbeta, Nkx, Nky, Onebody, Eri)
    modelparams2 = InputModel(nalpha + 1, nbeta - 1, Nkx, Nky, Onebody, Eri)
    elist = DiagonalizeAllMomentum(modelparams, 5)
    elist2 = DiagonalizeAllMomentum(modelparams2, 5)
    @show elist
    # plot the energy

    X = [[i for i in 1:Nkx*Nky] for j in 1:5]
    X = hcat(X...)

    plt.scatter(X, elist', s=40, alpha=1.0, c="blue", marker="x")
    plt.scatter(X, elist2', s=40, alpha=0.5, c="red", marker="o")
    plt.show()
    plt.close()


end

let
    nalpha = 3
    nbeta = 3
    Nkx = 6
    Nky = 6
    @time test(nalpha, nbeta, Nkx, Nky)

end
