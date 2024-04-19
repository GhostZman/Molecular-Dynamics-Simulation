//
//  MD Simulation.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 3/29/24.
//

import Foundation
import Observation
import SwiftUI

@Observable class MDSimulation: ObservableObject{
    var boxLength: Int = 1
    var numAtoms: Int = 4
    var nMax: Int = 4
    
    var atoms: [Particle2D] = []

    init() {
        self.nMax = self.numAtoms
        
    }
    
    //let objectWillChange = PassthroughSubject<Void, Never>()
    
    func runSimulation() {
        var atomIndex: Int
        var timeIndex: Int
        let nStep: Int = 5000
        let nPrint: Int = 100
        let nDim: Int = 1
        
        let stepSize: Double = 0.004
        var halfTimeStep: Double
        var PE: Double
        var KE: Double
        var temperature: Double
        let initialTemperature: Double = 10.0
        
        self.boxLength = Int(pow(Double(numAtoms), 1.0/Double(nDim)))
        self.numAtoms = Int(pow(Double(boxLength),Double(nDim)))
        
        print("numAtoms = \(numAtoms) boxLength = \(boxLength)")
        atomIndex = -1
        for ix in stride(from: 0, through: boxLength-1, by: 1){  //Set up lattice of side boxLength
            atomIndex = atomIndex+1
            let atomPosition = [Double(ix), 0.0]   // Inital velocities
            var atomVelocity = (Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0))/12.0 - 0.5
            var direction = Double.random(in: 0.0 ... Double.pi)
            atomVelocity = atomVelocity*sqrt(initialTemperature)   // Scale v with temperature
            var velocityVector = [atomVelocity*cos(direction), atomVelocity*sin(direction)]
            atoms.append(Particle2D(position: atomPosition, velocity: velocityVector,  boxSize: Double(boxLength)))
            
            print("init atomVelocity = \(atoms[atomIndex].velocity)")
        }
        halfTimeStep = stepSize/2.0
        timeIndex = 0
        KE = 0.0    // inital KE & PE
        PE = 0.0
        PE = updateForces()
        for atomIndex in stride(from: 0, through: numAtoms - 1, by: 1) {
            KE = KE + (pow(atoms[atomIndex].velocityMagnitude(), 2))/2
        }
        //print("\(timeIndex) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
        print("step: \(timeIndex)\nposition    velocity         force ")
        atoms.forEach { particle in
            print(particle.position, particle.velocity, particle.force)
        }
        
        for timeIndex in stride(from: 1, to: nStep, by: 1){     // Main loop
            usleep(5000)
            objectWillChange.send()
            for atomIndex in stride(from: 0, through: numAtoms - 1, by: 1) {   // Velocity Verlet
                for dim in stride(from: 0, through: 1, by: 1){
                    var newPosition = atoms[atomIndex].position[dim] + stepSize*(atoms[atomIndex].velocity[dim] + halfTimeStep*atoms[atomIndex].force[dim][1])      //PBC
                    if (newPosition <= 0.0){
                        newPosition = newPosition + Double(boxLength)
                    }
                    if (newPosition >= Double(boxLength)){
                        newPosition = newPosition - Double(boxLength)
                    }
                    atoms[atomIndex].changePosition(newPosition: newPosition, dimIndex: dim)
                }
            }
            PE = updateForces()
            KE = 0
            for atomIndex in stride(from: 0, through: numAtoms - 1, by: 1){
                for dim in stride(from: 0, through: 1, by: 1){
                    atoms[atomIndex].velocity[dim] = atoms[atomIndex].velocity[dim] + halfTimeStep*(atoms[atomIndex].force[dim][1] + atoms[atomIndex].force[dim][0])
                }
                KE = KE + (atoms[atomIndex].velocityMagnitude()*atoms[atomIndex].velocityMagnitude())/2
            }
            temperature = 2.0 * KE / (3.0 * Double(numAtoms))
            print("step: \(timeIndex)\nposition    velocity         force ")
            atoms.forEach { particle in
                print(particle.position, particle.velocity, particle.force)
            }
            if (timeIndex % nPrint == 0){
                //print("\(timeIndex) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
                /*print("step: \(timeIndex)\nposition    velocity         force ")
                atoms.forEach { particle in
                    print(particle.position, particle.velocity, particle.force)
                }
                */
            }
        }
                
    }
    
    func updateForces() -> Double {   //run this AFTER updating positions BEFORE updating velocities
        var forceOfiOnj: Double
        var separation2: Double
        var inverseSeparation2: Double = 0.0
        var displacementX: Double
        var displacementY: Double
        var angle: Double
        let cutoffSeparation2: Double = 9.0
        
        var newPE: Double = 0.0     // Initialize
        for i in stride(from: 0, through: numAtoms - 1, by: 1) {
            for dim in stride(from: 0, through: 1, by: 1) {
                atoms[i].force[dim][0] = atoms[i].force[dim][1]
                atoms[i].force[dim][1] = 0
            }
        }
        
        for i in stride(from: 0, through: numAtoms - 2, by: 1){
            for j in stride(from: i+1, through: numAtoms - 1, by: 1) {
                for atomImage in atoms[j].imagePositions{
                    displacementX = atoms[i].position[0] - atomImage[0]
                    displacementY = atoms[i].position[1] - atomImage[1]
                    if displacementX > 0 {
                        angle = atan(displacementY/displacementX)
                    }
                    else {
                        angle = atan(displacementY/displacementX) + Double.pi
                    }
                    separation2 = pow(displacementX, 2) + pow(displacementY, 2)
                    if (separation2 < cutoffSeparation2) {   // Cut off
                        if (abs(separation2) < 1.0E-7) {
                            separation2 = 1.0E-7
                        }
                        inverseSeparation2 = 1.0/separation2
                        forceOfiOnj = 48*(pow(inverseSeparation2, 3.0) - 0.5)*pow(inverseSeparation2, 3.0)
                        forceOfiOnj = forceOfiOnj*inverseSeparation2*sqrt(separation2)
                        
                        atoms[i].force[0][1] = atoms[i].force[0][1] + forceOfiOnj*cos(angle)
                        atoms[j].force[0][1] = atoms[j].force[0][1] - forceOfiOnj*cos(angle)
                        
                        atoms[i].force[1][1] = atoms[i].force[1][1] + forceOfiOnj*sin(angle)
                        atoms[j].force[1][1] = atoms[j].force[1][1] - forceOfiOnj*sin(angle)
                        
                        newPE = newPE + 4 * pow(inverseSeparation2, 3) * (pow(inverseSeparation2, 3) - 1)
                    }
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
