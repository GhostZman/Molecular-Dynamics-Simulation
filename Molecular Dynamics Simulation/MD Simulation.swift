//
//  MD Simulation.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 3/29/24.
//

import Foundation
import Observation

@Observable class MDSimulation{
    var L: Int = 1
    var numAtoms: Int = 8
    var nMax: Int = 513 //rename this once you figure out what it does
    var x: [Double] = []
    var fx: [[Double]] = [[]]

    init() {
        self.x = Array(repeating: 0, count: self.nMax)
        self.fx = Array(repeating: Array(repeating: 0, count: nMax), count: 2)
    }
    
    func forces(t: Int, PE: Double) -> Double {
        var fijx: Double
        var r2: Double
        var invr2: Double = 0.0
        var dx: Double
        var r2cut: Double = 9.0
        
        var newPE: Double = 0.0
        for i in stride(from: 0, to: numAtoms - 1, by: 1){
            fx[i][t] = 0
        }
        
        for i in stride(from: 0, through: numAtoms - 2, by: 1){
            for j in stride(from: i+1, through: numAtoms - 1, by: 1) {
                dx = x[i] - x[j]
                if (abs(dx) > 0.5*Double(L)){
                    dx = dx - sign(a: Double(L), b: dx)
                }
                r2 = dx * dx
                if (r2 < r2cut) {
                    if (r2 == 0.0) {
                        r2 = 0.0001
                    }
                    invr2 = 1.0/r2
                    fijx = 48*(pow(invr2, 3) - 0.5)*pow(invr2, 3)
                    fijx = fijx*invr2*dx
                    fx[i][t] = fx[i][t] + fijx
                    fx[j][t] = fx[j][t] - fijx
                    newPE = newPE + 4 * pow(invr2, 3) * (pow(invr2, 3) - 1)
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
