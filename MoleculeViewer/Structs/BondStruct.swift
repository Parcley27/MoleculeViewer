//
//  BondStruct.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 28/1/25.
//

import SwiftUI
import SceneKit
import Foundation

struct Bond {
    var atom1: Atom
    var atom2: Atom
    var type: Int // Bond type (e.g., single = 1, double = 2, etc.)
    var dipoleVector: SCNVector3? // Vector representing the dipole moment for this bond
    
}
