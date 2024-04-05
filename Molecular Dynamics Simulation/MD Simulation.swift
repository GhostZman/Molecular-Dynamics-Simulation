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
    var nAtom: Int = 4
    var nMax: Int = 4 //rename this once you figure out what it does
    var x: [Double] = []
    var fx: [[Double]] = [[]]

    init() {
        self.x = Array(repeating: 0, count: self.nMax)
        self.fx = Array(repeating: Array(repeating: 0, count: 2), count: nMax)
    }
    
    func runSimulation() {
        var t1: Int
        var t2: Int
        var i: Int
        var iTemp: Int
        var t: Int
        var nStep: Int = 5000
        var nPrint: Int = 100
        var nDim: Int = 1
        
        var h: Double = 0.0004
        var hover2: Double
        var PE: Double
        var KE: Double
        var T: Double
        var tInit: Double = 10.0
        var vx: [Double] = Array(repeating: 0, count: nMax)
        
        self.L = Int(pow(Double(nAtom), 1.0/Double(nDim)))
        self.nAtom = Int(pow(Double(L),Double(nDim)))
        
        print("nAtom = \(nAtom) L = \(L)")
        i = -1
        for ix in stride(from: 0, through: L-1, by: 1){  //Set up lattice of side L
            i = i+1
            x[i] = Double(ix)   // Inital velocities
            vx[i] = (Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0))/12.0 - 0.5
            vx[i] = vx[i]*sqrt(tInit)   // Scale v with temperature
            print("init vx = \(vx[i])")
        }
        t1 = 0      //t, t+h indicies
        t2 = 1
        hover2 = h/2.0
        t = 0
        KE = 0.0    // inital KE & PE
        PE = 0.0
        PE = forces(t: t1, PE: PE)
        for i in stride(from: 0, through: nAtom - 1, by: 1) {
            KE = KE + (vx[i]*vx[i])/2
        }
        //print("\(t) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
        print("\(t) \npos: \(x) \nvel: \(vx)")
        for t in stride(from: 1, to: nStep, by: 1){     // Main loop
            for i in stride(from: 0, through: nAtom - 1, by: 1) {   // Velocity Verlet
                PE = forces(t: t1, PE: PE)
                x[i] = x[i] + h*(vx[i] + hover2*fx[i][t1])      //PBC
                if (x[i] <= 0.0){
                    x[i] = x[i] + Double(L)
                }
                if (x[i] >= Double(L)){
                    x[i] = x[i] - Double(L)
                }
            }
            PE = forces(t: t2, PE: PE)
            KE = 0
            for i in stride(from: 0, through: nAtom - 1, by: 1){
                vx[i] = vx[i] + hover2*(fx[i][t1] + fx[i][t2])
                KE = KE + (vx[i]*vx[i])/2
            }
            T = 2.0 * KE / (3.0 * Double(nAtom))
            print("\(t) \npos: \(x) \nvel: \(vx) \nforces: \(fx)")
            if (t % nPrint == 0){
                //print("\(t) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
                //print("\(t) \npos: \(x) \nvel: \(vx)")
            }
            iTemp = t1      // Time t and t+h
            t1 = t2
            t2 = iTemp
        }
                
    }
    
    func forces(t: Int, PE: Double) -> Double {
        var fijx: Double
        var r2: Double
        var invr2: Double = 0.0
        var dx: Double
        var r2cut: Double = 9.0
        
        var newPE: Double = 0.0     // Initialize
        for i in stride(from: 0, to: nAtom - 1, by: 1){
            fx[i][t] = 0
        }
        
        for i in stride(from: 0, through: nAtom - 2, by: 1){
            for j in stride(from: i+1, through: nAtom - 1, by: 1) {
                dx = x[i] - x[j]
                if (abs(dx) > 0.5*Double(L)){   // PBC
                    dx = dx - sign(a: Double(L), b: dx)
                }
                r2 = dx * dx
                if (r2 < r2cut) {   // Cut off
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
