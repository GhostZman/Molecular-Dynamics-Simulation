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
    @ObservedObject var mySim = MDSimulation()
    @State private var points: [CGPoint] = []
    @State private var isDragging: Bool = false
    @State private var draggedIndex: Int?
    
    
    var body: some View {
        VStack {
            TimelineView(.animation) { timelineContext in
                let value = secondsValue(for: timelineContext.date)
                
                Canvas(
                    opaque: true,
                    colorMode: .linear,
                    rendersAsynchronously: false
                ) { context, size in
                    let newSize = size.applying(.init(scaleX: value, y:1))
                    let rect = CGRect(origin: .zero, size: newSize)
                    
                    context.fill(
                        Rectangle().path(in: rect),
                        with: .color(.red)
                    )
                }
            }
            Button("Start", action: sim)
        }
    }
    func sim(){
        mySim.runSimulation()
    }
    
    func secondsValue(for date: Date) -> Double{
        let seconds = Calendar.current.component(.second, from: date)
        return Double(seconds) / 60
    }
}

