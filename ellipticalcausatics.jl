@everywhere using Cubature
using Iterators
using PyCall
pygui(:gtk)
using PyPlot

# defining global variables
@everywhere a = 1.0
@everywhere b = 1.0
@everywhere c = 1.0

acclevel = 1
st = 10

k = zeros(Float64, 2)

@everywhere function sigma(x::Array{Float64,1})
    rho = 0.23873241463

    z = deepcopy(x)
    z[1] = z[1]/a
    z[2] = z[2]/b
    r = norm(z)

    if r<1.0
        2.0*c*rho*sqrt(1-r^2)
    else
        0.0
    end
end

@everywhere function dalpha(x::Array{Float64,1}, v::Array{Float64,1})
    m = ((k-x)/(norm(k-x))^2)*sigma(x)
    v[1] = m[1]
    v[2] = m[2]

    global k
    kd = k
    if isnan(v[1])
        k[1] += 1e-10
        dalpha(x,v)
    elseif isnan(isnan(v[2]))
        k[2] += 1e-10
        dalpha(x,v)
    end
    k = kd
end

@everywhere function quad(x::Array{Float64,1})
    global k
    k = x
    tuple((hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1])...)
    # hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1]
end

@everywhere function squad(x::Array{Float64,1})
    global k
    k = x
    tuple((k - hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1])...)
    # k - hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1]
end

@everywhere function squad(x::Float64, y::Float64)
    global k
    k = [x;y]
    tuple((k - hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1])...)
    # k - hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1]
end


function findCriticalCurves()
    acclevel = 2
    bnd = Any[Any[[-2.0;2.0],[-2.0;2.0]]]
    st = 16
    tbdry = cell(0)

    for acc=1:acclevel
        tbdry = cell(0)

        for bdry in bnd
            stepx = (bdry[1][2] - bdry[1][1])/(st-1)
			stepy = (bdry[2][2] - bdry[2][1])/(st-1)
            x = linspace(bdry[1][1]-stepx, bdry[1][2]+stepx, st+2)
            y = linspace(bdry[2][1]-stepy, bdry[2][2]+stepy, st+2)

            gd = Array{Float64,1}[]
            for ind in product(x,y)
                push!(gd, [ind...])
            end

            al = reshape(pmap(quad, gd), st+2, st+2)
            #@show al
            de = zeros(st, st)

            al11 = zeros(st, st)
            al12 = zeros(st, st)
            al21 = zeros(st, st)
            al22 = zeros(st, st)

            for i=2:st+1
                for j=2:st+1
                    al11[i-1,j-1] = (al[i+1,j][1] - al[i-1,j][1])/(2*stepx)
                    al12[i-1,j-1] = (al[i,j+1][1] - al[i,j-1][1])/(2*stepy)
                    al21[i-1,j-1] = (al[i+1,j][2] - al[i-1,j][2])/(2*stepx)
                    al22[i-1,j-1] = (al[i,j+1][2] - al[i,j-1][2])/(2*stepy)
                end
            end

            x = x[2:st+1]
            y = y[2:st+1]
            de = (1.0-al11).*(1.0-al22) - (((al12+al21)/2).^2)
            de = de./abs(de)
            # @show de

            for i=1:st-1
                for j=1:st-1
                    if de[i,j]*de[i+1,j]<1.0 || de[i,j]*de[i,j+1]<1.0 || de[i,j]*de[i+1,j+1]<1.0
                        push!(tbdry, Any[[x[i]-0.1*stepx;x[i+1]+0.1*stepx],[y[j]-0.1*stepy;y[j+1]+0.1*stepy]])
                    end
                end
            end
        end
        bnd = tbdry
    end

    fpointsx = [bdry[1][1] for bdry in bnd]
    fpointsy = [bdry[2][1] for bdry in bnd]
    
    scatter(fpointsx, fpointsy, color="red", marker="+")
    savefig("criticallines_2.svg")
    causatics = pmap(squad, fpointsx, fpointsy)
    cpointsx = [cau[1] for cau in causatics]
    cpointsy = [cau[2] for cau in causatics]
    scatter(cpointsx, cpointsy, color="blue", marker="o")
    savefig("causatics_2.svg")
end

@time findCriticalCurves()

function plotcircle()
    r = 1.0
    th = linspace(0.0, 2pi, 50)
    rx = r*cos(th)
    ry = r*sin(th)
    b = pmap(squad, rx, ry)
    bx = [bb[1] for bb in b]
    by = [bb[2] for bb in b]
    scatter(bx, by)
    savefig("testsquad.svg")
end

# @time plotcircle()
