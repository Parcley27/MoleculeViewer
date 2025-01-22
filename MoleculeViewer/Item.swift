//
//  Item.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import Foundation
import SwiftData

@Model
final class Item {
    var timestamp: Date
    
    init(timestamp: Date) {
        self.timestamp = timestamp
    }
}
