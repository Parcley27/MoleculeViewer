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
    func loadMolFile(named filename: String) -> String? {
        if let filePath = Bundle.main.path(forResource: filename, ofType: "mol") {
            return try? String(contentsOfFile: filePath, encoding: .utf8)
            
        }
        
        return nil
    }
    
    var body: some View {
        if let molFileContent = loadMolFile(named: "Sulfur Dioxide") {
            
            MoleculeView(molFile: molFileContent)
                .edgesIgnoringSafeArea(.all)
            
        }
    }
}

#Preview {
    ContentView()
        .modelContainer(for: Item.self, inMemory: true)
}
