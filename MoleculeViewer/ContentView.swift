//
//  ContentView.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import SwiftUI
import SwiftData

/*
struct ContentView: View {
    @Environment(\.modelContext) private var modelContext
    @Query private var items: [Item]

    var body: some View {
        NavigationSplitView {
            List {
                ForEach(items) { item in
                    NavigationLink {
                        Text("Item at \(item.timestamp, format: Date.FormatStyle(date: .numeric, time: .standard))")
                    } label: {
                        Text(item.timestamp, format: Date.FormatStyle(date: .numeric, time: .standard))
                    }
                }
                .onDelete(perform: deleteItems)
            }
            .toolbar {
                ToolbarItem(placement: .navigationBarTrailing) {
                    EditButton()
                }
                ToolbarItem {
                    Button(action: addItem) {
                        Label("Add Item", systemImage: "plus")
                    }
                }
            }
        } detail: {
            Text("Select an item")
        }
    }

    private func addItem() {
        withAnimation {
            let newItem = Item(timestamp: Date())
            modelContext.insert(newItem)
        }
    }

    private func deleteItems(offsets: IndexSet) {
        withAnimation {
            for index in offsets {
                modelContext.delete(items[index])
            }
        }
    }
}
 */

struct ContentView: View {
    var body: some View {
        ZStack {
            MoleculeView()
            
            VStack {
                Spacer()
                
                ZStack {
                    Capsule()
                        .fill(.white)
                        .opacity(0.5)
                        .blur(radius: 30)
                        .frame(width: 250, height: 100)
                    
                    VStack {
                        Text("C\u{2088}H\u{2081}\u{2080}N\u{2084}O\u{2082}")
                            .font(.system(size: 32, weight: .bold, design: .default))
                            .foregroundStyle(.black)
                        
                        Text("Caffine")
                            .font(.system(size: 24, weight: .bold, design: .default))
                            .foregroundStyle(.black)
                    }
                }
            }
        }
    }
}

#Preview {
    ContentView()
    
}
