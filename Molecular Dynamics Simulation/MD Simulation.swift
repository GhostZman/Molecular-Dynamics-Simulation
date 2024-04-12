//
//  MD Simulation.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 3/29/24.
//

import Foundation
import Observation

@Observable class MDSimulation{
    var boxLength: Int = 1
    var numAtoms: Int = 4
    var nMax: Int = 4
    var atomPostions: [Double] = []
    var atomForces: [[Double]] = [[]]
    
    var atoms: [Particle] = []

    init() {
        self.nMax = self.numAtoms
        self.atomPostions = Array(repeating: 0, count: self.nMax)
        self.atomForces = Array(repeating: Array(repeating: 0, count: 2), count: nMax)
        
    }
    
    func runSimulation() {
        var atomIndex: Int
        var timeIndex: Int
        var nStep: Int = 5000
        var nPrint: Int = 100
        var nDim: Int = 1
        
        var stepSize: Double = 0.004
        var halfTimeStep: Double
        var PE: Double
        var KE: Double
        var temperature: Double
        var initialTemperature: Double = 10.0
        var atomVelocities: [Double] = Array(repeating: 0, count: nMax)
        
        self.boxLength = Int(pow(Double(numAtoms), 1.0/Double(nDim)))
        self.numAtoms = Int(pow(Double(boxLength),Double(nDim)))
        
        print("numAtoms = \(numAtoms) boxLength = \(boxLength)")
        atomIndex = -1
        for ix in stride(from: 0, through: boxLength-1, by: 1){  //Set up lattice of side boxLength
            atomIndex = atomIndex+1
            let atomPostion = Double(ix)   // Inital velocities
            var atomVelocity = (Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0))/12.0 - 0.5
            atomVelocity = atomVelocity*sqrt(initialTemperature)   // Scale v with temperature
            
            atoms.append(Particle(position: atomPostion, velocity: atomVelocity, dimensions: nDim, boxSize: Double(boxLength)))
            
            print("init atomVelocity = \(atoms[atomIndex].velocity)")
        }
        halfTimeStep = stepSize/2.0
        timeIndex = 0
        KE = 0.0    // inital KE & PE
        PE = 0.0
        PE = updateForces()
        for atomIndex in stride(from: 0, through: numAtoms - 1, by: 1) {
            KE = KE + (pow(atoms[atomIndex].velocity, 2))/2
        }
        //print("\(timeIndex) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
        print("step: \(timeIndex)\nposition    velocity         force ")
        atoms.forEach { particle in
            print(particle.position, particle.velocity, particle.force)
        }
        
        for timeIndex in stride(from: 1, to: nStep, by: 1){     // Main loop
            for atomIndex in stride(from: 0, through: numAtoms - 1, by: 1) {   // Velocity Verlet
                atoms[atomIndex].position = atoms[atomIndex].position + stepSize*(atoms[atomIndex].velocity + halfTimeStep*atoms[atomIndex].force[1])      //PBC
                if (atoms[atomIndex].position <= 0.0){
                    atoms[atomIndex].position = atoms[atomIndex].position + Double(boxLength)
                }
                if (atoms[atomIndex].position >= Double(boxLength)){
                    atoms[atomIndex].position = atoms[atomIndex].position - Double(boxLength)
                }
            }
            PE = updateForces()
            KE = 0
            for atomIndex in stride(from: 0, through: numAtoms - 1, by: 1){
                atoms[atomIndex].velocity = atoms[atomIndex].velocity + halfTimeStep*(atoms[atomIndex].force[1] + atoms[atomIndex].force[0])
                KE = KE + (atoms[atomIndex].velocity*atoms[atomIndex].velocity)/2
            }
            temperature = 2.0 * KE / (3.0 * Double(numAtoms))
            print("step: \(timeIndex)\nposition    velocity         force ")
            atoms.forEach { particle in
                print(particle.position, particle.velocity, particle.force)
            }
            if (timeIndex % nPrint == 0){
                //print("\(timeIndex) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
                //print("\(timeIndex) \npos: \(atomPostions) \nvel: \(atomVelocities)")
            }
        }
                
    }
    
    func updateForces() -> Double {   //run this AFTER updating positions BEFORE updating velocities
        var forceOfiOnj: Double
        var separation2: Double
        var inverseSeparation2: Double = 0.0
        var displacement: Double
        let cutoffSeparation2: Double = 9.0
        
        var newPE: Double = 0.0     // Initialize
        for i in stride(from: 0, through: numAtoms - 1, by: 1){
            atoms[i].force[0] = atoms[i].force[1]
            atoms[i].force[1] = 0
        }
        
        for i in stride(from: 0, through: numAtoms - 2, by: 1){
            for j in stride(from: i+1, through: numAtoms - 1, by: 1) {
                displacement = atoms[i].position - atoms[j].position
                if (abs(displacement) > 0.5*Double(boxLength)){   // PBC
                    displacement = displacement - sign(a: Double(boxLength), b: displacement)
                }
                separation2 = displacement * displacement
                if (separation2 < cutoffSeparation2) {   // Cut off
                    if (abs(separation2) < 1.0E-7) {
                        separation2 = 1.0E-7
                    }
                    inverseSeparation2 = 1.0/separation2
                    forceOfiOnj = 48*(pow(inverseSeparation2, 3.0) - 0.5)*pow(inverseSeparation2, 3.0)
                    forceOfiOnj = forceOfiOnj*inverseSeparation2*displacement
                    
                    atoms[i].force[1] = atoms[i].force[1] + forceOfiOnj
                    atoms[j].force[1] = atoms[j].force[1] - forceOfiOnj
                    //print(atomForces[i][t],atomForces[j][t],forceOfiOnj,old1,old2)
                    newPE = newPE + 4 * pow(inverseSeparation2, 3) * (pow(inverseSeparation2, 3) - 1)
                }
            }
        }
        return newPE
    }
    func sign(a: Double, b: Double) -> Double{
        if (b >= 0) {
            return abs(a)
        }
        else {
            return -abs(a)
        }
    }
    
}
