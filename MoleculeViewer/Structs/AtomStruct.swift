//
//  AtomStruct.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 28/1/25.
//

import SwiftUI
import SceneKit
import Foundation

struct Atom {
    var id: Int
    var position: SCNVector3
    var element: String
    var size: CGFloat
    var color: UIColor
    var bonds: [Bond] // Bonds directly associated with the atom
    
}
