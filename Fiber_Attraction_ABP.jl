using Random, Distributions
using DelimitedFiles 
include("Construct_Fiber.jl")

#*************Set Parameters**************
const STEP = 10000  # Number of movement steps
const TICK = 10     #Number of steps between each save file
const NUMOFSIM = 5 #Number of simulation

const LENGTH = 1     #The length of box
const DENSITY = 0.3  # The density of system
const RADIUS = 0.02  #The radius of particles
const DELTATIME = 1  #particles[i].position[0] += DELTATIME * particles[i].velocity[0]
const V = 0.001     #Initial velocity

#const DIFFCOEFF = 1         #The diffusion coefficient
#const KECM = RADIUS * 0.0001 #The coefficient for ECM force
const ANGLE = 5                 #For the tolerance, velocities anti-parallel aligned within a prescribed tolerance

const PI = 3.1415926
const AREAOFBOX = LENGTH * LENGTH
const NUMBER = floor(Integer, AREAOFBOX * DENSITY / (PI * RADIUS * RADIUS))

#*************Add environment effect
const NUMGRID = 5     # Number of grids in one row/column
fiberChangeCoeff = 0  #for change angle of fiber. The coeff should be in the range [0,1]
const fiberRecoverCoeff = 0 # for change angle of fiber recovery. The coeff should be in the range [0,1]
const CHANGECOEFF = 0.25 #for each simulation, fiberRecoverCoeff += CHANGECOEFF
const fiberEffectCoeff = 0.5  #for change of velocity of particles. The coeff should be in the range [0,1]

#***********reconstruction of fiber
sizeBin = 0.5  # deltaR = sizeBin * r
rStepCoeff = 0.25 # rstep = initr * rStepCoeff
alpha = -20 # theoretical = exp.(alpha * correlationX)
const numOfRun = 10000 

#***********Random Parameters
normal = Normal(0, PI/64) #generate normal ditribution random number. (mean, std)
uniform = Uniform(0, 1) #generate uniform distribution random number 

mutable struct Grid
    direction :: Array{Float64,1}
    numOfParticles ::Int64
    particlesInGrid ::Array{Int64,1}
    position :: Array{Float64,1}
    initDirection :: Array{Float64,1}
end

mutable struct Particle
     position :: Array{Float64,1}
     velocity :: Array{Float64,1}
     labelofgrid :: Array{Int64,1}
     force :: Array{Float64,1}
     correctMove :: Array{Float64,1}
end

function boundaryDisplacement(a, b)
    delta = a - b
    if abs(delta) > LENGTH / 2
        if delta > 0
            delta -= LENGTH
        else
            delta += LENGTH
        end
    end
    return delta
end

function boundaryDistance(a, b)
    deltax = boundaryDisplacement(a.position[1], b.position[1])
    deltay = boundaryDisplacement(a.position[2], b.position[2])

    return sqrt(deltax * deltax + deltay * deltay)
end

function boundaryCondition(particle)
    if particle.position[1] >= LENGTH
        particle.position[1] = particle.position[1] - LENGTH
    end
    if particle.position[1] < 0
        particle.position[1] = particle.position[1] + LENGTH
    end
    if particle.position[2] >= LENGTH
        particle.position[2] = particle.position[2] - LENGTH
    end
    if particle.position[2] < 0
        particle.position[2] = particle.position[2] + LENGTH
    end
    return particle.position
end

function isOverlap(a, b)
    return (boundaryDistance(a, b) < 2 * RADIUS)
end

function isOppositeMove(particleA, particleB)
    d = boundaryDistance(particleA, particleB);
    dv = sqrt(particleA.velocity[1] * particleA.velocity[1] + particleA.velocity[2] * particleA.velocity[2])
    dv2 = sqrt(particleB.velocity[1] * particleB.velocity[1] + particleB.velocity[2] * particleB.velocity[2])

    deltax = boundaryDisplacement(particleA.position[1], particleB.position[1])
    deltay = boundaryDisplacement(particleA.position[2], particleB.position[2])
    deltax1 = boundaryDisplacement(particleB.position[1], particleA.position[1])
    deltay1 = boundaryDisplacement(particleB.position[2], particleA.position[2])

    cosTheta = (-deltax * particleA.velocity[1] + (-deltay) * particleA.velocity[2]) / (d * dv)
    cosBeta = (-deltax1 * particleB.velocity[1] + (-deltay1) * particleB.velocity[2]) / (d * dv2)

    return (cosTheta >= cos(ANGLE * PI / 180) && (cosBeta >= cos(ANGLE * PI / 180)) && (particleA.velocity[1] * particleB.velocity[1] + particleA.velocity[2] * particleB.velocity[2] < 0))
end

function initializationGrids()
    grids = Array{Grid}(undef, NUMGRID,NUMGRID)
    Threads.@threads for i = 1:NUMGRID
        for j = 1:NUMGRID
            randDir = rand(uniform)
            #cos(2*PI*randDir),sin(2*PI*randDir)
            #grids[i,j] = Grid([1,0],0,[],[(i-1)*(LENGTH / NUMGRID) + LENGTH / NUMGRID / 2, (j-1)*(LENGTH / NUMGRID) + LENGTH / NUMGRID / 2],[1,0])
            grids[i,j] = Grid([cos(2*PI*randDir),sin(2*PI*randDir)],0,[],[(i-1)*(LENGTH / NUMGRID) + LENGTH / NUMGRID / 2, (j-1)*(LENGTH / NUMGRID) + LENGTH / NUMGRID / 2],[cos(2*PI*randDir),sin(2*PI*randDir)])
        end
    end

    return grids
end

function initializationParticles()
    particles= Array{Particle}(undef, NUMBER)
    Threads.@threads for ii = 1:NUMBER
        particles[ii] = Particle([0,0],[0,0],[0,0],[0,0],[0,0])
    end
    i = 1
    while i <= NUMBER
        particles[i].position = [LENGTH * rand(uniform) , LENGTH * rand(uniform)] 
        particles[i].labelofgrid = [min(ceil(Integer, particles[i].position[1] / (LENGTH / NUMGRID)), NUMGRID) , min(ceil(Integer, particles[i].position[2] / (LENGTH / NUMGRID)), NUMGRID)]
  
        hasOverlap = false
        for j = 1:i
            if i != j && isOverlap(particles[i], particles[j])
                hasOverlap = true
                break
            end
        end
        if !hasOverlap
            i += 1
        end
    end
#pragma omp parallel for private(ii)
    Threads.@threads for ii = 1:NUMBER
        particles[ii].velocity = [V * LENGTH * (rand(uniform)-0.5) , V * LENGTH * (rand(uniform)-0.5)]
    end
    return particles
end

function rotateAngleAtoB(vectorA::Array, vectorB::Array, coeff::Number)
    vectorAangle = atan(vectorA[2], vectorA[1])
    vectorBangle = atan(vectorB[2], vectorB[1])
    diff = vectorBangle - vectorAangle
    sign = (diff < 0) ? -1.0 : 1.0
    while abs(diff) > PI / 2
        diff -= PI * sign
    end
    return diff * coeff
end

function rotateAtoB(vectorA::Array, vectorB::Array, coeff::Number)
    absVectorA = sqrt(vectorA[1] * vectorA[1] + vectorA[2] * vectorA[2])
    vectorAangle = atan(vectorA[2], vectorA[1])
    vectorAangle += rotateAngleAtoB(vectorA, vectorB, coeff)
    vectorA = [absVectorA * cos(vectorAangle),absVectorA * sin(vectorAangle)]
    return vectorA
end

function rotateAtoB(vectorA::Array, rotateAngle::Number)
    absVectorA = sqrt(vectorA[1] * vectorA[1] + vectorA[2] * vectorA[2])
    vectorAangle = atan(vectorA[2], vectorA[1])
    vectorAangle += rotateAngle 
    vectorA = [absVectorA * cos(vectorAangle),absVectorA * sin(vectorAangle)]
    return vectorA
end


function fiberAffected(particles, grids)
#pragma omp parallel for private(i)

    #Threads.@threads 
    for i = 1: NUMBER
        particles[i].labelofgrid = [min(ceil(Integer, particles[i].position[1] / (LENGTH / NUMGRID)), NUMGRID) , min(ceil(Integer, particles[i].position[2] / (LENGTH / NUMGRID)), NUMGRID)]

        grids[particles[i].labelofgrid[1], particles[i].labelofgrid[2]].numOfParticles +=1 
        append!(grids[particles[i].labelofgrid[1], particles[i].labelofgrid[2]].particlesInGrid, i)  
    end

    #pragma omp parallel for private(ii)
    #Threads.@threads 
    for i = 1:NUMGRID
        for j = 1:NUMGRID
            rotateAngle = 0
            if (grids[i, j].numOfParticles >= 2)
                for p1 = 1:grids[i, j].numOfParticles
                     for p2 = p1 + 1 :  grids[i, j].numOfParticles
                         a = grids[i, j].particlesInGrid[p1]
                         b = grids[i, j].particlesInGrid[p2]
                         if (isOppositeMove(particles[a], particles[b]))
                            rotateAngle += rotateAngleAtoB(grids[i, j].direction, particles[a].position - particles[b].position, fiberChangeCoeff)
                         end
                    end
                end
                grids[i, j].direction = rotateAtoB(grids[i, j].direction, rotateAngle)    
            else
                grids[i, j].direction = rotateAtoB(grids[i, j].direction, grids[i, j].initDirection, fiberRecoverCoeff);
            end
            grids[i, j].numOfParticles = 0
        end
    end
end

function moveParticles(particles, grids)
    # Move all particles. Update the positions of particles 
    # according to the current velocities.Velocities are updated 
    # by Gaussian random walk. Gaussian distribution is set by
    # global normal distribution generator.

    #pragma omp parallel for private(i)
    Threads.@threads for i = 1:NUMBER 
        particles[i].labelofgrid = [min(ceil(Integer, particles[i].position[1] / (LENGTH / NUMGRID)), NUMGRID) , min(ceil(Integer, particles[i].position[2] / (LENGTH / NUMGRID)), NUMGRID)]
        #Add force
        particles[i].velocity[1] += particles[i].force[1]
        particles[i].velocity[2] += particles[i].force[2]
        #environment effect
        particles[i].velocity = rotateAtoB(particles[i].velocity, grids[particles[i].labelofgrid[1], particles[i].labelofgrid[2]].direction, fiberEffectCoeff)
        #rotation diffusion
        particles[i].velocity = rotateAtoB(particles[i].velocity, rand(normal))
        #Update positions
        particles[i].position[1] += DELTATIME * particles[i].velocity[1];
        particles[i].position[2] += DELTATIME * particles[i].velocity[2];
        # Correct the position according to periodic boundary condition.
        particles[i].position = boundaryCondition(particles[i])
    end
end

function correctOverlap(particleA, particleB)
    
    # Offset overlap between two particles.If two particles overlap, 
    # move each one-half the overlap distance along their center-to-center axis.
    
    d = boundaryDistance(particleA, particleB)
    deltax = boundaryDisplacement(particleA.position[1], particleB.position[1])
    deltay = boundaryDisplacement(particleA.position[2], particleB.position[2])
    deltaD = 2 * RADIUS - d

    particleA.correctMove[1] += 0.5 * deltaD * deltax / d
    particleA.correctMove[2] += 0.5 * deltaD * deltay / d

end

function recordCorrectOverlap(particles)
#pragma omp parallel for private(j)
    Threads.@threads for i = 1:NUMBER
        for j = 1: NUMBER
            if j != i
                if (isOverlap(particles[i], particles[j]))
                    correctOverlap(particles[i], particles[j])
                end
            end
        end
    end
end

function moveCorrectOverlap(particles)
#pragma omp parallel for private(i)
    recordCorrectOverlap(particles)
    for i = 1:NUMBER
        particles[i].position[1] += particles[i].correctMove[1]
        particles[i].position[2] += particles[i].correctMove[2]
        particles[i].correctMove = [0, 0]
        particles[i].position = boundaryCondition(particles[i])
    end
end


function simulation(simID)
    # Run a complete simulation from start. simID labels the save detination for data.
    # After initialization, particles move for a number of steps specified by the global
    # constant STEP. Motion data is saved for every TICK steps.
    grids = initializationGrids()
    generateFiber(grids,sizeBin, rStepCoeff,alpha,numOfRun)
    particles = initializationParticles()
    
    positiondata = Array{Float64}(undef, NUMBER, 4)
    fiberPositiondata = Array{Float64}(undef, NUMGRID*NUMGRID,2)
    for i = 1:NUMBER
        positiondata[i,:] = [particles[i].position[1], particles[i].position[2], particles[i].velocity[1], particles[i].velocity[2]]
    end
    writedlm("ChangeCoeff$(round(fiberChangeCoeff; digits=2))/position0.csv", positiondata, ',')
    k = 1
    for i = 1:NUMGRID
        for j = 1:NUMGRID
            fiberPositiondata[k,:] =[grids[i,j].direction[1],grids[i,j].direction[2]]
            k +=1 
        end
    end
    writedlm("ChangeCoeff$(round(fiberChangeCoeff; digits=2))/fiberPosition0.csv",  fiberPositiondata, ',')

    for step = 1:STEP
        for i = 1:NUMBER
            particles[i].force=[0,0]
        end
        moveParticles(particles, grids)
        moveCorrectOverlap(particles)
        fiberAffected(particles, grids)

        if step % TICK == 0
            for i = 1:NUMBER
                positiondata[i,:] = [particles[i].position[1], particles[i].position[2], particles[i].velocity[1], particles[i].velocity[2]]
            end
            writedlm("ChangeCoeff$(round(fiberChangeCoeff; digits=2))/position$(step).csv", positiondata, ',')

            k = 1
            for i = 1:NUMGRID
                for j = 1:NUMGRID
                    fiberPositiondata[k,:] =[grids[i,j].direction[1],grids[i,j].direction[2]]
                    k +=1 
                end
            end
            writedlm("ChangeCoeff$(round(fiberChangeCoeff; digits=2))/fiberPosition$(step).csv",  fiberPositiondata, ',')
        end
    end 
end

function runSimulation()
    for i = 1:NUMOFSIM
        println("simulation: ", i)
        simulation(i)
        global fiberChangeCoeff += CHANGECOEFF
    end
end

@time runSimulation()


