//
//  ContentView.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 3/29/24.
//

//import SwiftUI
//
//struct ContentView: View {
//    let mySim = MDSimulation()
//    var body: some View {
//        VStack {
//            Canvas { context, size in
//                context.stroke(Path(roundedRect: CGRect(origin: .zero, size: size), cornerSize: CGSize(width: 20, height: 50)), with: .color(.red))
//                context.stroke(Path(roundedRect: CGRect(origin: .zero, size: size), cornerSize: CGSize(width: 50, height: 20)), with: .color(.green))
//            }
//            .frame(width: 300, height: 200)
//            
//            Button("Start", action: sim)
//            
//        }
//        .background(.black)
//        .ignoresSafeArea()
//        .padding()
//    }
//    
//    func sim(){
//        mySim.runSimulation()
//    }
//}
import SwiftUI

struct ContentView: View {
    let mySim = MDSimulation()
    @State private var points: [CGPoint] = []
    @State private var isDragging: Bool = false
    @State private var draggedIndex: Int?
    
    var body: some View {
        VStack {
            DrawingView(atoms: mySim.atoms)
                .frame(maxWidth: 100, maxHeight: 100)
            
            Button("Start", action: sim)
        }
    }
    func sim(){
        mySim.runSimulation()
    }
}

