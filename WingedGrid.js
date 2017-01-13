/* WARNING, CURRENTLY CAN ONLY SUBDIVIDE TRIANGULAR GRIDS. SUBDIVIDING DUALS IS UNDEFINED */
var WingedEdge = (function () {
    function WingedEdge(sourceGrid, edgeIndex) {
        this.grid = sourceGrid;
        this.index = edgeIndex;
    }
    WingedEdge.prototype.firstVertexA = function () {
        return this.grid.edges_firstVertexA[this.index];
    };
    WingedEdge.prototype.firstVertexB = function () {
        return this.grid.edges_firstVertexB[this.index];
    };
    WingedEdge.prototype.faceA = function () {
        return this.grid.edges_faceA[this.index];
    };
    WingedEdge.prototype.faceB = function () {
        return this.grid.edges_faceB[this.index];
    };
    WingedEdge.prototype.prevA = function () {
        return this.grid.edges_prevA[this.index];
    };
    WingedEdge.prototype.nextA = function () {
        return this.grid.edges_nextA[this.index];
    };
    WingedEdge.prototype.prevB = function () {
        return this.grid.edges_prevB[this.index];
    };
    WingedEdge.prototype.nextB = function () {
        return this.grid.edges_nextB[this.index];
    };
    /********* Base **********/
    WingedEdge.prototype.nextEdgeForFace = function (faceIndex) {
        if (this.faceA() == faceIndex) {
            return [this.nextA(), null];
        }
        else if (this.faceB() == faceIndex) {
            return [this.nextB(), null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.prevEdgeForFace = function (faceIndex) {
        if (this.faceA() == faceIndex) {
            return [this.prevA(), null];
        }
        else if (this.faceB() == faceIndex) {
            return [this.prevB(), null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.firstVertexForFace = function (faceIndex) {
        if (this.faceA() == faceIndex) {
            return [this.firstVertexA(), null];
        }
        else if (this.faceB() == faceIndex) {
            return [this.firstVertexB(), null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.secondVertexForFace = function (faceIndex) {
        if (this.faceA() == faceIndex) {
            return [this.firstVertexB(), null];
        }
        else if (this.faceB() == faceIndex) {
            return [this.firstVertexA(), null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.nextEdgeForVertex = function (vertexIndex) {
        if (this.firstVertexA() == vertexIndex) {
            return this.prevEdgeForFace(this.faceA());
        }
        else if (this.firstVertexB() == vertexIndex) {
            return this.prevEdgeForFace(this.faceB());
        }
        return [-1, Error("Edge not associated with vertex.")];
    };
    WingedEdge.prototype.prevEdgeForVertex = function (vertexIndex) {
        if (this.firstVertexA() == vertexIndex) {
            return this.nextEdgeForFace(this.faceB());
        }
        else if (this.firstVertexB() == vertexIndex) {
            return this.nextEdgeForFace(this.faceA());
        }
        return [-1, Error("Edge not associated with vertex.")];
    };
    WingedEdge.prototype.adjacentForFace = function (faceIndex) {
        if (this.faceA() == faceIndex) {
            return [this.faceB(), null];
        }
        else if (this.faceB() == faceIndex) {
            return [this.faceA(), null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.adjacentForVertex = function (vertexIndex) {
        if (this.firstVertexA() == vertexIndex) {
            return [this.firstVertexB(), null];
        }
        else if (this.firstVertexB() == vertexIndex) {
            return [this.firstVertexA(), null];
        }
        return [-1, Error("Edge not associated with vertex.")];
    };
    return WingedEdge;
}());
var WingedFace = (function () {
    function WingedFace(sourceGrid, faceIndex) {
        this.grid = sourceGrid;
        this.index = faceIndex;
    }
    WingedFace.prototype.edges = function () {
        var end;
        if (this.index == this.grid.faces_edgeOffsets.length - 1) {
            end = this.grid.faces_edges.length;
        }
        else {
            end = this.grid.faces_edgeOffsets[this.index + 1];
        }
        var edges = new Array(end - this.grid.faces_edgeOffsets[this.index]);
        var i = 0;
        for (var faceEdge = this.grid.faces_edgeOffsets[this.index]; faceEdge < end; ++faceEdge) {
            edges[i] = new WingedEdge(this.grid, this.grid.faces_edges[faceEdge]);
            i++;
        }
        return edges;
    };
    return WingedFace;
}());
var WingedVertex = (function () {
    function WingedVertex(sourceGrid, faceIndex) {
        this.grid = sourceGrid;
        this.index = faceIndex;
    }
    WingedVertex.prototype.coords = function () {
        return [this.grid.vertices_coords[this.index * 3], this.grid.vertices_coords[this.index * 3 + 1], this.grid.vertices_coords[this.index * 3 + 2]];
    };
    ;
    WingedVertex.prototype.edges = function () {
        var end;
        if (this.index == this.grid.vertices_edgeOffsets.length - 1) {
            end = this.grid.vertices_edges.length;
        }
        else {
            end = this.grid.vertices_edgeOffsets[this.index + 1];
        }
        var edges = new Array(end - this.grid.vertices_edgeOffsets[this.index]);
        var i = 0;
        for (var faceEdge = this.grid.vertices_edgeOffsets[this.index]; faceEdge < end; ++faceEdge) {
            edges[i] = new WingedEdge(this.grid, this.grid.vertices_edges[faceEdge]);
            i++;
        }
        return edges;
    };
    WingedVertex.prototype.vertexNeighbors = function () {
        var neighborsIndecies = this.grid.neighborsForVertex(this.index);
        var neighbors = new Array(neighborsIndecies.length);
        for (var i = neighborsIndecies.length - 1; i >= 0; i--) {
            neighbors[i] = new WingedVertex(this.grid, neighborsIndecies[i]);
        }
        return neighbors;
    };
    ;
    return WingedVertex;
}());
var WingedGrid = (function () {
    function WingedGrid(faceCount, edgeCount, vertCount, totalFaceEdgeCount, totalVertexEdgeCount) {
        this.faces_edgeOffsets = new Int32Array(faceCount);
        this.faces_edges = new Int32Array(totalFaceEdgeCount);
        this.edges_firstVertexA = new Int32Array(edgeCount);
        this.edges_firstVertexA = new Int32Array(edgeCount);
        this.edges_firstVertexB = new Int32Array(edgeCount);
        this.edges_faceA = new Int32Array(edgeCount);
        this.edges_faceB = new Int32Array(edgeCount);
        this.edges_prevA = new Int32Array(edgeCount);
        this.edges_nextA = new Int32Array(edgeCount);
        this.edges_prevB = new Int32Array(edgeCount);
        this.edges_nextB = new Int32Array(edgeCount);
        this.vertices_coords = new Float32Array(vertCount * 3);
        this.vertices_edges = new Int32Array(totalVertexEdgeCount);
        this.vertices_edgeOffsets = new Int32Array(vertCount);
    }
    WingedGrid.prototype.faces = function (index) {
        return new WingedFace(this, index);
    };
    WingedGrid.prototype.edges = function (index) {
        return new WingedEdge(this, index);
    };
    WingedGrid.prototype.vertices = function (index) {
        return new WingedVertex(this, index);
    };
    /********* Base **********/
    // TODO, OPTIMIZE FOR NEW STRUCTURE
    WingedGrid.prototype.neighborsForFace = function (faceIndex) {
        var neighbors = Array();
        for (var i = this.faces(faceIndex).edges().length - 1; i >= 0; i--) {
            var _a = this.faces(faceIndex).edges()[i].adjacentForFace(faceIndex), neighbor = _a[0], err = _a[1];
            neighbors[i] = neighbor;
        }
        return neighbors;
    };
    WingedGrid.prototype.neighborsForVertex = function (vertexIndex) {
        var neighbors = Array();
        for (var i = this.vertices(vertexIndex).edges().length - 1; i >= 0; i--) {
            var _a = this.vertices(vertexIndex).edges()[i].adjacentForVertex(vertexIndex), neighbor = _a[0], error = _a[1];
            neighbors[i] = neighbor;
        }
        return neighbors;
    };
    /*********** SPHERE ************/
    WingedGrid.prototype.normalizeVerticesToDistanceFromOrigin = function (wantedLength) {
        for (var i = this.vertices_edgeOffsets.length - 1; i >= 0; i--) {
            var vertex = this.vertices(i);
            var currentLength = vectorLength(vertex.coords());
            this.vertices_coords[i * 3 + 0] = vertex.coords()[0] * wantedLength / currentLength;
            this.vertices_coords[i * 3 + 1] = vertex.coords()[1] * wantedLength / currentLength;
            this.vertices_coords[i * 3 + 2] = vertex.coords()[2] * wantedLength / currentLength;
        }
    };
    /******** DUAL *********/
    WingedGrid.prototype.createDual = function () {
        var dualGrid = new WingedGrid(this.vertices_edgeOffsets.length, this.edges_firstVertexA.length, this.faces_edgeOffsets.length, 0, 0); // swap vert and face counts
        // create faces, just swap with old vert edges
        dualGrid.faces_edges = this.vertices_edges;
        dualGrid.faces_edgeOffsets = this.vertices_edgeOffsets;
        // create vertices
        // swap with old face edges
        dualGrid.vertices_edges = this.faces_edges;
        dualGrid.vertices_edgeOffsets = this.faces_edgeOffsets;
        for (var i = this.faces_edgeOffsets.length - 1; i >= 0; i--) {
            var faceEdges = this.faces(i).edges();
            // set coords from center of old face
            var faceCenter = [0, 0, 0];
            var count = 0;
            for (var j = faceEdges.length - 1; j >= 0; j--) {
                var edge = faceEdges[j];
                var _a = edge.firstVertexForFace(i), vertexIndex = _a[0], err = _a[1];
                if (err != null) {
                    return [dualGrid, err];
                }
                var vertex = this.vertices(vertexIndex);
                faceCenter[0] += vertex.coords()[0];
                faceCenter[1] += vertex.coords()[1];
                faceCenter[2] += vertex.coords()[2];
                count += 1;
            }
            dualGrid.vertices_coords[i * 3 + 0] = faceCenter[0] / count;
            dualGrid.vertices_coords[i * 3 + 1] = faceCenter[1] / count;
            dualGrid.vertices_coords[i * 3 + 2] = faceCenter[2] / count;
        }
        // set edges, swap with old
        dualGrid.edges_faceA = this.edges_firstVertexA;
        dualGrid.edges_faceB = this.edges_firstVertexB;
        dualGrid.edges_firstVertexA = this.edges_faceB;
        dualGrid.edges_firstVertexB = this.edges_faceA;
        // set prev and next
        for (var faceIndex = dualGrid.faces.length - 1; faceIndex >= 0; faceIndex--) {
            var faceEdges = dualGrid.faces(faceIndex).edges();
            for (var faceEdgeIndex = faceEdges.length - 1; faceEdgeIndex >= 0; faceEdgeIndex--) {
                var edgeIndex = faceEdges[faceEdgeIndex].index;
                var theEdge = dualGrid.edges(edgeIndex);
                if (theEdge.faceA() == faceIndex) {
                    if (faceEdgeIndex == 0) {
                        dualGrid.edges_prevA[edgeIndex] = faceEdges[faceEdges.length - 1].index;
                        dualGrid.edges_nextA[edgeIndex] = faceEdges[1].index;
                    }
                    else if (faceEdgeIndex == faceEdges.length - 1) {
                        dualGrid.edges_prevA[edgeIndex] = faceEdges[faceEdges.length - 1].index;
                        dualGrid.edges_nextA[edgeIndex] = faceEdges[0].index;
                    }
                    else {
                        dualGrid.edges_prevA[edgeIndex] = faceEdges[faceEdgeIndex - 1].index;
                        dualGrid.edges_nextA[edgeIndex] = faceEdges[faceEdgeIndex + 1].index;
                    }
                }
                if (theEdge.faceB() == faceIndex) {
                    if (faceEdgeIndex == 0) {
                        dualGrid.edges_prevB[edgeIndex] = faceEdges[faceEdges.length - 1].index;
                        dualGrid.edges_nextB[edgeIndex] = faceEdges[1].index;
                    }
                    else if (faceEdgeIndex == faceEdges.length - 1) {
                        dualGrid.edges_prevB[edgeIndex] = faceEdges[faceEdges.length - 1].index;
                        dualGrid.edges_nextB[edgeIndex] = faceEdges[0].index;
                    }
                    else {
                        dualGrid.edges_prevB[edgeIndex] = faceEdges[faceEdgeIndex - 1].index;
                        dualGrid.edges_nextB[edgeIndex] = faceEdges[faceEdgeIndex + 1].index;
                    }
                }
            }
        }
        return [dualGrid, null];
    };
    /************* Subdivision helpers ****************/
    // rename vertexInexAtClockwiseIndexOnFaceToSubdivide
    WingedGrid.prototype.vertexIndexAtClockwiseIndexOnOldFace = function (faceIndex, edgeInFaceIndex, clockwiseVertexIndex, edgeSubdivisions) {
        var edge = this.faces(faceIndex).edges()[edgeInFaceIndex];
        var edgeIndex = edge.index;
        if (edge.faceA() == faceIndex) {
            return this.vertices_edgeOffsets.length + edgeIndex * edgeSubdivisions + clockwiseVertexIndex;
        }
        if (edge.faceB() == faceIndex) {
            return this.vertices_edgeOffsets.length + edgeIndex * edgeSubdivisions + edgeSubdivisions - 1 - clockwiseVertexIndex;
        }
        return -1;
    };
    WingedGrid.prototype.edgeIndexAtClockwiseIndexOnOldFace = function (faceIndex, edgeInFaceIndex, clockwiseEdgeIndex, edgeSubdivisions) {
        var oldEdge = this.faces(faceIndex).edges()[edgeInFaceIndex];
        var oldEdgeIndex = oldEdge.index;
        if (oldEdge.faceA() == faceIndex) {
            return oldEdgeIndex * (edgeSubdivisions + 1) + clockwiseEdgeIndex;
        }
        if (oldEdge.faceB() == faceIndex) {
            return oldEdgeIndex * (edgeSubdivisions + 1) + edgeSubdivisions - clockwiseEdgeIndex;
        }
        return -1;
    };
    /*************** SUBDIVISIONS ***************/
    WingedGrid.prototype.subdivideTriangles = function (edgeSubdivisions) {
        if (edgeSubdivisions < 1) {
            return [null, Error("Invalid number of subdivisions.")];
        }
        // subdividing each edge n times produces a number of faces equal
        //  to the base face count multiplied by (1/2(n+2)(n+1) + 1/2(n+1)(n))
        var faceCount = this.faces_edgeOffsets.length * ((edgeSubdivisions + 2) * (edgeSubdivisions + 1) / 2 + (edgeSubdivisions + 1) * edgeSubdivisions / 2);
        var dividedGrid = new WingedGrid(faceCount, 
        // since each face 'owns' 1/2 of three edges, there are 1.5 times as
        //  many edges as faces
        3 * faceCount / 2, 
        // Euler-somebody or other gives us the vertex count of
        //  1/2*(face count) + 2
        faceCount / 2 + 2, faceCount * 3, // three for each face
        (faceCount / 2 + 2 - this.vertices_edgeOffsets.length) * 6 + this.vertices_edges.length);
        // Invalidate all values that will be set later, useful for error checking
        for (var i = dividedGrid.faces_edges.length - 1; i >= 0; i--) {
            dividedGrid.faces_edges[i] = -1;
        }
        for (var i = dividedGrid.faces_edgeOffsets.length - 1; i >= 0; i--) {
            dividedGrid.faces_edgeOffsets[i] = -1;
        }
        for (var i = dividedGrid.edges_firstVertexA.length - 1; i >= 0; i--) {
            dividedGrid.edges_faceA[i] = -1;
            dividedGrid.edges_faceB[i] = -1;
            dividedGrid.edges_firstVertexA[i] = -1;
            dividedGrid.edges_firstVertexB[i] = -1;
            dividedGrid.edges_nextA[i] = -1;
            dividedGrid.edges_nextB[i] = -1;
            dividedGrid.edges_prevA[i] = -1;
            dividedGrid.edges_prevB[i] = -1;
        }
        // vertices corrisponding with old ones will have the same number of
        //  associated edges to preserve the Euler characteristic ( =2 for S2)
        for (var i = dividedGrid.vertices_edges.length - 1; i >= 0; i--) {
            dividedGrid.vertices_edges[i] = -1;
        }
        var currentOffset = 0;
        for (var i = 0; i < dividedGrid.vertices_edgeOffsets.length; ++i) {
            if (i < this.vertices_edgeOffsets.length) {
                currentOffset = this.vertices_edgeOffsets[i];
            }
            else {
                currentOffset += 6;
            }
            dividedGrid.vertices_edgeOffsets[i] = currentOffset;
        }
        for (var i = dividedGrid.vertices_coords.length - 1; i >= 0; i--) {
            dividedGrid.vertices_coords[i] = Infinity;
        }
        /***************** Subdivide the grid ****************/
        /******************* EDGE SUBDIVISION ***********************/
        /****** Old Edge Subdivision goes in the first section of the new array, ordered by edge *******/
        var origVertexCount = this.vertices_edgeOffsets.length;
        var origEdgeCount = this.edges_faceA.length;
        for (var i_1 = origEdgeCount - 1; i_1 >= 0; i_1--) {
            var edge_1 = this.edges(i_1);
            // first edge has origional vertex
            dividedGrid.edges_firstVertexA[i_1 * (edgeSubdivisions + 1)] = edge_1.firstVertexA();
            dividedGrid.edges_firstVertexB[i_1 * (edgeSubdivisions + 1)] = origVertexCount + i_1 * edgeSubdivisions;
            for (var j = 1; j < edgeSubdivisions; j++) {
                // set the edge vertex indecies
                dividedGrid.edges_firstVertexA[i_1 * (edgeSubdivisions + 1) + j] = origVertexCount + i_1 * edgeSubdivisions + j - 1;
                dividedGrid.edges_firstVertexB[i_1 * (edgeSubdivisions + 1) + j] = origVertexCount + i_1 * edgeSubdivisions + j;
            }
            // connect last new edge
            dividedGrid.edges_firstVertexA[i_1 * (edgeSubdivisions + 1) + edgeSubdivisions] = origVertexCount + i_1 * edgeSubdivisions + edgeSubdivisions - 1;
            dividedGrid.edges_firstVertexB[i_1 * (edgeSubdivisions + 1) + edgeSubdivisions] = edge_1.firstVertexB();
        }
        /********* Edges created interior to old faces go in the second section, ordered by face. ****/
        for (var faceIndex_1 = this.faces_edgeOffsets.length - 1; faceIndex_1 >= 0; faceIndex_1--) {
            var oldFace = this.faces(faceIndex_1);
            // vertex offset for vertices interior to the face
            var vertexOffset_1 = origVertexCount + origEdgeCount * edgeSubdivisions + (edgeSubdivisions * (edgeSubdivisions - 1) / 2) * faceIndex_1;
            var edgeOffset = (edgeSubdivisions + 1) * origEdgeCount + 3 * (edgeSubdivisions) * (edgeSubdivisions + 1) / 2 * faceIndex_1;
            if (edgeSubdivisions == 1) {
                dividedGrid.edges_firstVertexA[edgeOffset + 0] = origVertexCount + oldFace.edges()[0].index;
                dividedGrid.edges_firstVertexB[edgeOffset + 0] = origVertexCount + oldFace.edges()[1].index;
                dividedGrid.edges_firstVertexA[edgeOffset + 1] = origVertexCount + oldFace.edges()[0].index;
                dividedGrid.edges_firstVertexB[edgeOffset + 1] = origVertexCount + oldFace.edges()[2].index;
                dividedGrid.edges_firstVertexA[edgeOffset + 2] = origVertexCount + oldFace.edges()[1].index;
                dividedGrid.edges_firstVertexB[edgeOffset + 2] = origVertexCount + oldFace.edges()[2].index;
            }
            else {
                // first row
                dividedGrid.edges_firstVertexA[edgeOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1, edgeSubdivisions);
                dividedGrid.edges_firstVertexB[edgeOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, 0, edgeSubdivisions);
                dividedGrid.edges_firstVertexA[edgeOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1, edgeSubdivisions);
                dividedGrid.edges_firstVertexB[edgeOffset + 1] = vertexOffset_1 + 0;
                dividedGrid.edges_firstVertexA[edgeOffset + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, 0, edgeSubdivisions);
                dividedGrid.edges_firstVertexB[edgeOffset + 2] = vertexOffset_1 + 0;
                // middle rows
                var rowOffset;
                var i = 0;
                for (i = 1; i < edgeSubdivisions - 1; i++) {
                    rowOffset = i * (i + 1) * 3 / 2;
                    // first border
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1 - i, edgeSubdivisions);
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 0] = vertexOffset_1 + (i * (i - 1) / 2);
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1 - i, edgeSubdivisions);
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 1] = vertexOffset_1 + (i * (i - 1) / 2) + i;
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 2] = vertexOffset_1 + (i * (i - 1) / 2);
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 2] = vertexOffset_1 + (i * (i - 1)) / 2 + i;
                    // interior of face
                    for (var j = 1; j < i; j++) {
                        dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                        dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset_1 + (i * (i - 1) / 2) + j;
                        dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 1] = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                        dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 1] = vertexOffset_1 + (i * (i - 1) / 2) + i + j;
                        dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 2] = vertexOffset_1 + (i * (i - 1) / 2) + j;
                        dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 2] = vertexOffset_1 + (i * (i - 1) / 2) + i + j;
                    }
                    // second border
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 0] = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 1] = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 1] = vertexOffset_1 + (i * (i - 1) / 2) + 2 * i;
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 2] = vertexOffset_1 + (i * (i - 1) / 2) + 2 * i;
                }
                // last row
                i = edgeSubdivisions - 1; // should already be set, but just incase
                rowOffset = i * (i + 1) * 3 / 2;
                // border 1
                dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, 0, edgeSubdivisions);
                dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 0] = vertexOffset_1 + (i * (i - 1) / 2);
                dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, 0, edgeSubdivisions);
                dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - 1, edgeSubdivisions);
                dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 2] = vertexOffset_1 + (i * (i - 1) / 2);
                dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - 1, edgeSubdivisions);
                // middle
                for (var j = 1; j < i; j++) {
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset_1 + (i * (i - 1) / 2) + j;
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 1] = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - j - 1, edgeSubdivisions);
                    dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 2] = vertexOffset_1 + (i * (i - 1) / 2) + j;
                    dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - j - 1, edgeSubdivisions);
                }
                // border 2
                dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 0] = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 1] = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, 0, edgeSubdivisions);
                dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, 0, edgeSubdivisions);
            }
        }
        // create the faces, update the edges
        /******************* FACE SUBDIVISION ***********************/
        for (var faceIndex_2 = this.faces_edgeOffsets.length - 1; faceIndex_2 >= 0; faceIndex_2--) {
            // get the number of faces we divide this one into
            var subFaceCount = (edgeSubdivisions + 2) * (edgeSubdivisions + 1) / 2 + (edgeSubdivisions + 1) * edgeSubdivisions / 2;
            // number of internal edges created on this face
            //  equal to ( (edgeSubdivisions + 1) choose 2 ) * 3
            var subEdgeCount = 3 * edgeSubdivisions * (edgeSubdivisions + 1) / 2;
            // each new face set starts at
            var indexStart = subFaceCount * faceIndex_2;
            var edgeOffset = (edgeSubdivisions + 1) * origEdgeCount + subEdgeCount * faceIndex_2;
            // first corner
            // set edges in face
            dividedGrid.faces_edgeOffsets[indexStart] = (indexStart) * 3; // three edges per face
            dividedGrid.faces_edges[indexStart * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 0, edgeSubdivisions, edgeSubdivisions);
            dividedGrid.faces_edges[indexStart * 3 + 1] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 1, 0, edgeSubdivisions);
            dividedGrid.faces_edges[indexStart * 3 + 2] = edgeOffset;
            // loop through the middle section of faces
            // edges grow by 1/2*i*(i-1)*3
            for (var i_2 = 1; i_2 < edgeSubdivisions; i_2++) {
                var rowIndexStart = i_2 * (i_2 + 1) / 2 + i_2 * (i_2 - 1) / 2;
                // edge
                dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart] = (indexStart + rowIndexStart) * 3;
                dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 0, edgeSubdivisions - i_2, edgeSubdivisions);
                dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 1] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + 1;
                dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 2] = edgeOffset + i_2 * (i_2 + 1) * 3 / 2;
                // middle
                // up to one less than next index start
                for (var j = 1; j < (i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - rowIndexStart - 1; j++) {
                    if (j % 2 == 1) {
                        dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
                        dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + (j - 1) * 3 / 2;
                        dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + (j - 1) * 3 / 2 + 2;
                        dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + (j - 1) * 3 / 2 + 1;
                    }
                    else {
                        dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
                        dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = edgeOffset + i_2 * (i_2 + 1) * 3 / 2 + j * 3 / 2;
                        dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + (j - 2) * 3 / 2 + 2;
                        dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + j * 3 / 2 + 1;
                    }
                }
                // edge
                dividedGrid.faces_edgeOffsets[indexStart + (i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - 1] = (indexStart + (i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - 1) * 3;
                dividedGrid.faces_edges[(indexStart + (i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - 1) * 3 + 0] = edgeOffset + i_2 * (i_2 + 1) * 3 / 2 + ((i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - rowIndexStart - 1) * 3 / 2;
                dividedGrid.faces_edges[(indexStart + (i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - 1) * 3 + 1] = edgeOffset + i_2 * (i_2 - 1) * 3 / 2 + ((i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - rowIndexStart - 3) * 3 / 2 + 2;
                dividedGrid.faces_edges[(indexStart + (i_2 + 1) * (i_2 + 2) / 2 + (i_2 + 1) * i_2 / 2 - 1) * 3 + 2] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 1, i_2, edgeSubdivisions);
            }
            var rowIndexStart = edgeSubdivisions * (edgeSubdivisions + 1) / 2 + edgeSubdivisions * (edgeSubdivisions - 1) / 2;
            // bottom corner 1
            dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart] = (indexStart + rowIndexStart) * 3;
            dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 0, 0, edgeSubdivisions);
            dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + 1;
            dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 2] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 2, edgeSubdivisions, edgeSubdivisions);
            // bottom edge
            for (var j = 1; j < subFaceCount - rowIndexStart - 1; j++) {
                if (j % 2 == 1) {
                    dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
                    dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2;
                    dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2 + 2;
                    dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2 + 1;
                }
                else {
                    dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
                    dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 2, edgeSubdivisions - j / 2, edgeSubdivisions);
                    dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 2) * 3 / 2 + 2;
                    dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + j * 3 / 2 + 1;
                }
            }
            // bottom corner 2
            dividedGrid.faces_edgeOffsets[indexStart + subFaceCount - 1] = (indexStart + subFaceCount - 1) * 3;
            dividedGrid.faces_edges[(indexStart + subFaceCount - 1) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 2, 0, edgeSubdivisions);
            dividedGrid.faces_edges[(indexStart + subFaceCount - 1) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (subFaceCount - rowIndexStart - 3) * 3 / 2 + 2;
            dividedGrid.faces_edges[(indexStart + subFaceCount - 1) * 3 + 2] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 1, edgeSubdivisions, edgeSubdivisions);
        }
        // set edge faces from the previously build edge arrays
        for (var faceIndex = dividedGrid.faces_edgeOffsets.length - 1; faceIndex >= 0; faceIndex--) {
            var edges = dividedGrid.faces(faceIndex).edges();
            var edgesLength = edges.length;
            var thisEdge, nextEdge, prevEdge;
            var prevEdge_1 = edges[edgesLength - 1];
            var thisEdge_1 = edges[0];
            var nextEdge_1 = edges[1];
            // test vertices
            if (thisEdge_1.firstVertexA() == prevEdge_1.firstVertexA() || thisEdge_1.firstVertexA() == prevEdge_1.firstVertexB()) {
                // check the next edge also matches and face A is not set
                if (thisEdge_1.firstVertexB() == nextEdge_1.firstVertexA() || thisEdge_1.firstVertexB() == nextEdge_1.firstVertexB()) {
                    if (thisEdge_1.faceA() == -1) {
                        dividedGrid.edges_faceA[edges[0].index] = faceIndex;
                        dividedGrid.edges_prevA[edges[0].index] = edges[edgesLength - 1].index;
                        dividedGrid.edges_nextA[edges[0].index] = edges[1].index;
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge_1.index + " With edge set: " + edges);
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                }
            }
            else if (thisEdge_1.firstVertexB() == prevEdge_1.firstVertexA() || thisEdge_1.firstVertexB() == prevEdge_1.firstVertexB()) {
                // check the next edge also matches and face B is not set
                if (thisEdge_1.firstVertexA() == nextEdge_1.firstVertexA() || thisEdge_1.firstVertexA() == nextEdge_1.firstVertexB()) {
                    if (thisEdge_1.faceB() == -1) {
                        dividedGrid.edges_faceB[edges[0].index] = faceIndex;
                        dividedGrid.edges_prevB[edges[0].index] = edges[edgesLength - 1].index;
                        dividedGrid.edges_nextB[edges[0].index] = edges[1].index;
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge_1.index + " With edge set: " + edges);
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                }
            }
            else {
                console.log("For face " + faceIndex + ". Edges Don't share a vertex!");
            }
            // loop through middle edges
            for (var i_3 = 1; i_3 < edgesLength - 1; i_3++) {
                prevEdge_1 = edges[i_3 - 1];
                thisEdge_1 = edges[i_3];
                nextEdge_1 = edges[i_3 + 1];
                // test vertcies
                if (thisEdge_1.firstVertexA() == prevEdge_1.firstVertexA() || thisEdge_1.firstVertexA() == prevEdge_1.firstVertexB()) {
                    // check the next edge also matches
                    if (thisEdge_1.firstVertexB() == nextEdge_1.firstVertexA() || thisEdge_1.firstVertexB() == nextEdge_1.firstVertexB()) {
                        if (thisEdge_1.faceA() == -1) {
                            dividedGrid.edges_faceA[edges[i_3].index] = faceIndex;
                            dividedGrid.edges_prevA[edges[i_3].index] = edges[i_3 - 1].index;
                            dividedGrid.edges_nextA[edges[i_3].index] = edges[i_3 + 1].index;
                        }
                        else {
                            console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge_1.index + " With edge set: " + edges);
                        }
                    }
                    else {
                        console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                    }
                }
                else if (thisEdge_1.firstVertexB() == prevEdge_1.firstVertexA() || thisEdge_1.firstVertexB() == prevEdge_1.firstVertexB()) {
                    // check the next edge also matches
                    if (thisEdge_1.firstVertexA() == nextEdge_1.firstVertexA() || thisEdge_1.firstVertexA() == nextEdge_1.firstVertexB()) {
                        if (thisEdge_1.faceB() == -1) {
                            dividedGrid.edges_faceB[edges[i_3].index] = faceIndex;
                            dividedGrid.edges_prevB[edges[i_3].index] = edges[i_3 - 1].index;
                            dividedGrid.edges_nextB[edges[i_3].index] = edges[i_3 + 1].index;
                        }
                        else {
                            console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge_1.index + " With edge set: " + edges);
                        }
                    }
                    else {
                        console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Edges Don't share a vertex!");
                }
            }
            // last edge
            prevEdge_1 = edges[edgesLength - 2];
            thisEdge_1 = edges[edgesLength - 1];
            nextEdge_1 = edges[0];
            // test vertices
            if (thisEdge_1.firstVertexA() == prevEdge_1.firstVertexA() || thisEdge_1.firstVertexA() == prevEdge_1.firstVertexB()) {
                // check the next edge also matches
                if (thisEdge_1.firstVertexB() == nextEdge_1.firstVertexA() || thisEdge_1.firstVertexB() == nextEdge_1.firstVertexB()) {
                    if (thisEdge_1.faceA() == -1) {
                        dividedGrid.edges_faceA[edges[edgesLength - 1].index] = faceIndex;
                        dividedGrid.edges_prevA[edges[edgesLength - 1].index] = edges[edgesLength - 2].index;
                        dividedGrid.edges_nextA[edges[edgesLength - 1].index] = edges[0].index;
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge_1.index + " With edge set: " + edges);
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                }
            }
            else if (thisEdge_1.firstVertexB() == prevEdge_1.firstVertexA() || thisEdge_1.firstVertexB() == prevEdge_1.firstVertexB()) {
                // check the next edge also matches
                if (thisEdge_1.firstVertexA() == nextEdge_1.firstVertexA() || thisEdge_1.firstVertexA() == nextEdge_1.firstVertexB()) {
                    if (thisEdge_1.faceB() == -1) {
                        dividedGrid.edges_faceB[edges[edgesLength - 1].index] = faceIndex;
                        dividedGrid.edges_prevB[edges[edgesLength - 1].index] = edges[edgesLength - 2].index;
                        dividedGrid.edges_nextB[edges[edgesLength - 1].index] = edges[0].index;
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge_1.index + " With edge set: " + edges);
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                }
            }
            else {
                console.log("For face " + faceIndex + ". Edges Don't share a vertex!");
            } //*/
        }
        // create the vertices
        /******************* VERTEX SUBDIVISION ***********************/
        // set coords for the origional verts
        for (var i_4 = this.vertices_edgeOffsets.length - 1; i_4 >= 0; i_4--) {
            var vertex = this.vertices(i_4);
            dividedGrid.vertices_coords[i_4 * 3 + 0] = vertex.coords()[0];
            dividedGrid.vertices_coords[i_4 * 3 + 1] = vertex.coords()[1];
            dividedGrid.vertices_coords[i_4 * 3 + 2] = vertex.coords()[2];
        }
        // subdivide along each edge
        for (var i_5 = this.edges_faceA.length - 1; i_5 >= 0; i_5--) {
            var edge_2 = this.edges(i_5);
            var firstVertex_1 = dividedGrid.vertices(edge_2.firstVertexA());
            var secondVertex_1 = dividedGrid.vertices(edge_2.firstVertexB());
            // angle to subdivide
            var angleToSubdivide_1 = void 0;
            angleToSubdivide_1 = vectorAngle(firstVertex_1.coords(), secondVertex_1.coords());
            // angle between origin, first vertex, and second vertex
            var vectorA_1 = [-1, -1, -1], vectorB_1 = [-1, -1, -1];
            vectorA_1[0] = -1 * firstVertex_1.coords()[0];
            vectorA_1[1] = -1 * firstVertex_1.coords()[1];
            vectorA_1[2] = -1 * firstVertex_1.coords()[2];
            vectorB_1[0] = secondVertex_1.coords()[0] - firstVertex_1.coords()[0];
            vectorB_1[1] = secondVertex_1.coords()[1] - firstVertex_1.coords()[1];
            vectorB_1[2] = secondVertex_1.coords()[2] - firstVertex_1.coords()[2];
            var cornerAngle_1 = void 0;
            cornerAngle_1 = vectorAngle(vectorA_1, vectorB_1);
            // unit vector from first to second vertex
            var stepDirection_1 = [-1, -1, -1];
            stepDirection_1[0] = vectorB_1[0] / vectorLength(vectorB_1);
            stepDirection_1[1] = vectorB_1[1] / vectorLength(vectorB_1);
            stepDirection_1[2] = vectorB_1[2] / vectorLength(vectorB_1);
            // origional radius of the
            var sphereRadius = vectorLength(firstVertex_1.coords());
            var divisionLength_1 = void 0;
            for (var j = 0; j < edgeSubdivisions; j++) {
                // find the new vertex position and create the vertex
                // but don't correct it's length yet
                divisionLength_1 = Math.sin(angleToSubdivide_1 * ((j + 1) / (edgeSubdivisions + 1))) * sphereRadius / Math.sin(Math.PI - cornerAngle_1 - angleToSubdivide_1 * ((j + 1) / (edgeSubdivisions + 1)));
                dividedGrid.vertices_coords[(origVertexCount + i_5 * edgeSubdivisions + j) * 3 + 0] = firstVertex_1.coords()[0] + stepDirection_1[0] * divisionLength_1;
                dividedGrid.vertices_coords[(origVertexCount + i_5 * edgeSubdivisions + j) * 3 + 1] = firstVertex_1.coords()[1] + stepDirection_1[1] * divisionLength_1;
                dividedGrid.vertices_coords[(origVertexCount + i_5 * edgeSubdivisions + j) * 3 + 2] = firstVertex_1.coords()[2] + stepDirection_1[2] * divisionLength_1;
            }
        }
        // subdivide face interior
        for (var faceIndex_3 = 0; faceIndex_3 < this.faces_edgeOffsets.length; faceIndex_3++) {
            var vertexOffset = origVertexCount + origEdgeCount * edgeSubdivisions + edgeSubdivisions * (edgeSubdivisions - 1) / 2 * faceIndex_3;
            // only if we have more than one division
            if (edgeSubdivisions > 1) {
                for (var i_6 = 0; i_6 < edgeSubdivisions - 1; i_6++) {
                    var firstVertex = dividedGrid.vertices(this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_3, 0, edgeSubdivisions - 2 - i_6, edgeSubdivisions));
                    var secondVertex = dividedGrid.vertices(this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_3, 1, 1 + i_6, edgeSubdivisions));
                    // angle to subdivide
                    var angleToSubdivide;
                    angleToSubdivide = vectorAngle(firstVertex.coords(), secondVertex.coords());
                    // angle between origin, first vertex, and second vertex
                    var vectorA = [-1, -1, -1], vectorB = [-1, -1, -1];
                    vectorA[0] = -1 * firstVertex.coords()[0];
                    vectorA[1] = -1 * firstVertex.coords()[1];
                    vectorA[2] = -1 * firstVertex.coords()[2];
                    vectorB[0] = secondVertex.coords()[0] - firstVertex.coords()[0];
                    vectorB[1] = secondVertex.coords()[1] - firstVertex.coords()[1];
                    vectorB[2] = secondVertex.coords()[2] - firstVertex.coords()[2];
                    var cornerAngle;
                    cornerAngle = vectorAngle(vectorA, vectorB);
                    // unit vector from first to second vertex
                    var stepDirection = [-1, -1, -1];
                    stepDirection[0] = vectorB[0] / vectorLength(vectorB);
                    stepDirection[1] = vectorB[1] / vectorLength(vectorB);
                    stepDirection[2] = vectorB[2] / vectorLength(vectorB);
                    var firstVectorLength = vectorLength(firstVertex.coords());
                    var divisionLength = 0;
                    for (var j = 0; j < i_6 + 1; j++) {
                        divisionLength = Math.sin(angleToSubdivide * ((j + 1) / (i_6 + 2))) * firstVectorLength / Math.sin(Math.PI - cornerAngle - angleToSubdivide * ((j + 1) / (i_6 + 2)));
                        if (vertexOffset + (i_6 * (i_6 + 1) / 2) + j >= dividedGrid.vertices_edgeOffsets.length) {
                            console.log("breakpoint");
                        }
                        dividedGrid.vertices_coords[(vertexOffset + (i_6 * (i_6 + 1) / 2) + j) * 3 + 0] = firstVertex.coords()[0] + stepDirection[0] * divisionLength;
                        dividedGrid.vertices_coords[(vertexOffset + (i_6 * (i_6 + 1) / 2) + j) * 3 + 1] = firstVertex.coords()[1] + stepDirection[1] * divisionLength;
                        dividedGrid.vertices_coords[(vertexOffset + (i_6 * (i_6 + 1) / 2) + j) * 3 + 2] = firstVertex.coords()[2] + stepDirection[2] * divisionLength;
                    }
                }
            }
        }
        // set vertex edge array
        // loop through edges so we only have to touch each one once
        for (var index = dividedGrid.edges_faceA.length - 1; index >= 0; index--) {
            var edge = dividedGrid.edges(index);
            var nextEdgeIndex = -1;
            var nextEdge;
            var theVertex = dividedGrid.vertices(edge.firstVertexA());
            var theVertexIndex = theVertex.index;
            var error;
            if (dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[0]] == -1) {
                _a = edge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _a[0], error = _a[1];
                nextEdge = dividedGrid.edges(nextEdgeIndex);
                dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertex.index] + 0] = index;
                var i = 1;
                while (index != nextEdgeIndex) {
                    dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertex.index] + i] = nextEdgeIndex;
                    _b = nextEdge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _b[0], error = _b[1];
                    nextEdge = dividedGrid.edges(nextEdgeIndex);
                    i = i + 1;
                }
            }
            theVertex = dividedGrid.vertices(edge.firstVertexB());
            theVertexIndex = theVertex.index;
            if (dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[0]] == -1) {
                _c = edge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _c[0], error = _c[1];
                nextEdge = dividedGrid.edges(nextEdgeIndex);
                dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertex.index] + 0] = index;
                var i = 1;
                while (index != nextEdgeIndex) {
                    dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertex.index] + i] = nextEdgeIndex;
                    _d = nextEdge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _d[0], error = _d[1];
                    nextEdge = dividedGrid.edges(nextEdgeIndex);
                    i = i + 1;
                }
            }
        }
        return [dividedGrid, null];
        var _a, _b, _c, _d;
    };
    return WingedGrid;
}());
function baseIcosahedron() {
    var icosahedron = new WingedGrid(0, 0, 0, 0, 0); // arrays will be set manually
    icosahedron.faces_edges = Int32Array.from([
        // cap 1
        9, 4, 0 // 0
        ,
        1, 5, 0 // 1
        ,
        2, 6, 1 // 2
        ,
        3, 7, 2 // 3
        ,
        4, 8, 3 // 4
        ,
        14, 19, 10 // 5
        ,
        15, 11, 10 // 6
        ,
        16, 12, 11 // 7
        ,
        17, 13, 12 // 8
        ,
        18, 14, 13 // 9
        ,
        29, 28, 9 // 10
        ,
        29, 20, 17 // 11
        ,
        21, 20, 5 // 12
        ,
        21, 22, 18 // 13
        ,
        23, 22, 6 // 14
        ,
        23, 24, 19 // 15
        ,
        25, 24, 7 // 16
        ,
        25, 26, 15 // 17
        ,
        27, 26, 8 // 18
        ,
        27, 28, 16 // 19
    ]);
    // set offsets, three apiece
    icosahedron.faces_edgeOffsets = new Int32Array(20);
    for (var i = 0; i < 20; ++i) {
        icosahedron.faces_edgeOffsets[i] = i * 3;
    }
    var goldenRatio = 1.61803398875;
    icosahedron.vertices_coords = Float32Array.from([
        0, 1, goldenRatio,
        0, 1, -goldenRatio,
        0, -1, goldenRatio,
        0, -1, -goldenRatio,
        -goldenRatio, 0, 1,
        goldenRatio, 0, 1,
        -goldenRatio, 0, -1,
        goldenRatio, 0, -1,
        -1, goldenRatio, 0,
        -1, -goldenRatio, 0,
        1, goldenRatio, 0,
        1, -goldenRatio, 0
    ]);
    // Edges are duplicate info, could be done pragmatically after constucting
    // edges
    // Must be in counter-clockwise order
    // y-z plane rectangle
    icosahedron.vertices_edges = Int32Array.from([
        4, 3, 2, 1, 0,
        19, 24, 25, 15, 10,
        5, 20, 29, 9, 0,
        11, 12, 13, 14, 10,
        // x-z plane rectangle
        6, 22, 21, 5, 1,
        9, 28, 27, 8, 4,
        18, 22, 23, 19, 14,
        15, 26, 27, 16, 11,
        // x-y plane rectangle
        7, 24, 23, 6, 2,
        17, 20, 21, 18, 13,
        8, 26, 25, 7, 3,
        16, 28, 29, 17, 12
    ]);
    // set vertex offsets
    icosahedron.vertices_edgeOffsets = new Int32Array(12);
    for (var i = 0; i < 12; ++i) {
        icosahedron.vertices_edgeOffsets[i] = i * 5;
    }
    icosahedron.edges_firstVertexA = Int32Array.from([0, 0, 0, 0, 0, 4, 8, 10, 5, 2, 1, 7, 11, 9, 6, 1, 7, 11, 9, 6, 2, 4, 4, 8, 8, 10, 10, 5, 5, 2]);
    icosahedron.edges_firstVertexB = Int32Array.from([2, 4, 8, 10, 5, 2, 4, 8, 10, 5, 3, 3, 3, 3, 3, 7, 11, 9, 6, 1, 9, 9, 6, 6, 1, 1, 7, 7, 11, 11]);
    icosahedron.edges_faceA = Int32Array.from([0, 1, 2, 3, 4, 1, 2, 3, 4, 0, 5, 6, 7, 8, 9, 6, 7, 8, 9, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 10]);
    icosahedron.edges_faceB = Int32Array.from([1, 2, 3, 4, 0, 12, 14, 16, 18, 10, 6, 7, 8, 9, 5, 17, 19, 11, 13, 15, 12, 13, 14, 15, 16, 17, 18, 19, 10, 11]);
    icosahedron.edges_prevA = Int32Array.from([4, 0, 1, 2, 3, 1, 2, 3, 4, 0, 19, 15, 16, 17, 18, 10, 11, 12, 13, 14, 29, 5, 21, 6, 23, 7, 25, 8, 27, 9]);
    icosahedron.edges_nextA = Int32Array.from([9, 5, 6, 7, 8, 0, 1, 2, 3, 4, 14, 10, 11, 12, 13, 11, 12, 13, 14, 10, 17, 20, 18, 22, 19, 24, 15, 26, 16, 28]);
    icosahedron.edges_prevB = Int32Array.from([5, 6, 7, 8, 9, 20, 22, 24, 26, 28, 11, 12, 13, 14, 10, 26, 28, 20, 22, 24, 21, 18, 23, 19, 25, 15, 27, 16, 29, 17]);
    icosahedron.edges_nextB = Int32Array.from([1, 2, 3, 4, 0, 21, 23, 25, 27, 29, 15, 16, 17, 18, 19, 25, 27, 29, 21, 23, 5, 22, 6, 24, 7, 26, 8, 28, 9, 20]);
    return [icosahedron, null];
}
function vectorAngle(first, second) {
    return Math.acos((first[0] * second[0] + first[1] * second[1] + first[2] * second[2]) / (Math.sqrt(first[0] * first[0] + first[1] * first[1] + first[2] * first[2]) * Math.sqrt(second[0] * second[0] + second[1] * second[1] + second[2] * second[2])));
}
function vectorLength(vector) {
    return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}
