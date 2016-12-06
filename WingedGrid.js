var WingedEdge = (function () {
    function WingedEdge() {
    }
    /********* Base **********/
    WingedEdge.prototype.nextEdgeForFace = function (faceIndex) {
        if (this.faceA == faceIndex) {
            return [this.nextA, null];
        }
        else if (this.faceB == faceIndex) {
            return [this.nextB, null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.prevEdgeForFace = function (faceIndex) {
        if (this.faceA == faceIndex) {
            return [this.prevA, null];
        }
        else if (this.faceB == faceIndex) {
            return [this.prevB, null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.firstVertexForFace = function (faceIndex) {
        if (this.faceA == faceIndex) {
            return [this.firstVertexA, null];
        }
        else if (this.faceB == faceIndex) {
            return [this.firstVertexB, null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.secondVertexForFace = function (faceIndex) {
        if (this.faceA == faceIndex) {
            return [this.firstVertexB, null];
        }
        else if (this.faceB == faceIndex) {
            return [this.firstVertexA, null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.nextEdgeForVertex = function (vertexIndex) {
        if (this.firstVertexA == vertexIndex) {
            return this.prevEdgeForFace(this.faceA);
        }
        else if (this.firstVertexB == vertexIndex) {
            return this.prevEdgeForFace(this.faceB);
        }
        return [-1, Error("Edge not associated with vertex.")];
    };
    WingedEdge.prototype.prevEdgeForVertex = function (vertexIndex) {
        if (this.firstVertexA == vertexIndex) {
            return this.nextEdgeForFace(this.faceB);
        }
        else if (this.firstVertexB == vertexIndex) {
            return this.nextEdgeForFace(this.faceA);
        }
        return [-1, Error("Edge not associated with vertex.")];
    };
    WingedEdge.prototype.adjacentForFace = function (faceIndex) {
        if (this.faceA == faceIndex) {
            return [this.faceB, null];
        }
        else if (this.faceB == faceIndex) {
            return [this.faceA, null];
        }
        return [-1, Error("Edge not associated with face.")];
    };
    WingedEdge.prototype.adjacentForVertex = function (vertexIndex) {
        if (this.firstVertexA == vertexIndex) {
            return [this.firstVertexB, null];
        }
        else if (this.firstVertexB == vertexIndex) {
            return [this.firstVertexA, null];
        }
        return [-1, Error("Edge not associated with vertex.")];
    };
    return WingedEdge;
}());
var WingedFace = (function () {
    function WingedFace() {
        this.edges = new Array();
    }
    return WingedFace;
}());
var WingedVertex = (function () {
    function WingedVertex() {
        this.coords = [0, 0, 0];
        this.edges = new Array();
        this.vertexNeighbors = new Array();
    }
    return WingedVertex;
}());
var WingedGrid = (function () {
    function WingedGrid(faceCount, edgeCount, vertCount) {
        this.faces = new Array(faceCount);
        for (var i = this.faces.length - 1; i >= 0; i--) {
            this.faces[i] = new WingedFace();
        }
        this.edges = new Array(edgeCount);
        for (var i = this.edges.length - 1; i >= 0; i--) {
            this.edges[i] = new WingedEdge();
        }
        this.vertices = new Array(vertCount);
        for (var i = this.vertices.length - 1; i >= 0; i--) {
            this.vertices[i] = new WingedVertex();
        }
    }
    /********* Base **********/
    WingedGrid.prototype.neighborsForFace = function (faceIndex) {
        var neighbors = Array();
        for (var i = this.faces[faceIndex].edges.length - 1; i >= 0; i--) {
            var _a = this.edges[this.faces[faceIndex].edges[i]].adjacentForFace(faceIndex), neighbor = _a[0], err = _a[1];
            neighbors[i] = neighbor;
        }
        return neighbors;
    };
    WingedGrid.prototype.neighborsForVertex = function (vertexIndex) {
        var neighbors = Array();
        if (this.vertices[vertexIndex].vertexNeighbors.length == 0) {
            for (var i = this.vertices[vertexIndex].edges.length - 1; i >= 0; i--) {
                var _a = this.edges[this.vertices[vertexIndex].edges[i]].adjacentForVertex(vertexIndex), neighbor = _a[0], error = _a[1];
                this.vertices[vertexIndex].vertexNeighbors[i] = neighbor;
            }
        }
        return neighbors;
    };
    /*********** SPHERE ************/
    WingedGrid.prototype.normalizeVerticesToDistanceFromOrigin = function (wantedLength) {
        for (var i = this.vertices.length - 1; i >= 0; i--) {
            var vertex = this.vertices[i];
            var currentLength = vectorLength(vertex.coords);
            this.vertices[i].coords[0] = vertex.coords[0] * wantedLength / currentLength;
            this.vertices[i].coords[1] = vertex.coords[1] * wantedLength / currentLength;
            this.vertices[i].coords[2] = vertex.coords[2] * wantedLength / currentLength;
        }
    };
    /******** DUAL *********/
    WingedGrid.prototype.createDual = function () {
        var dualGrid = new WingedGrid(this.vertices.length, this.edges.length, this.faces.length); // swap vert and face counts
        // create faces
        for (var i = this.vertices.length - 1; i >= 0; i--) {
            var vertex = this.vertices[i];
            for (var j = vertex.edges.length - 1; j >= 0; j--) {
                dualGrid.faces[i].edges[j] = vertex.edges[j];
            }
        }
        // create vertices
        for (var i = this.faces.length - 1; i >= 0; i--) {
            var face = this.faces[i];
            // set edges
            for (var j = face.edges.length - 1; j >= 0; j--) {
                dualGrid.vertices[i].edges[j] = face.edges[j];
            }
            // set coords from center of old face
            var faceCenter = [0, 0, 0];
            var count = 0;
            for (var j = face.edges.length - 1; j >= 0; j--) {
                var edge = this.edges[face.edges[j]];
                var _a = edge.firstVertexForFace(i), vertexIndex = _a[0], err = _a[1];
                if (err != null) {
                    return [dualGrid, err];
                }
                var vertex = this.vertices[vertexIndex];
                faceCenter[0] += vertex.coords[0];
                faceCenter[1] += vertex.coords[1];
                faceCenter[2] += vertex.coords[2];
                count += 1;
            }
            dualGrid.vertices[i].coords[0] = faceCenter[0] / count;
            dualGrid.vertices[i].coords[1] = faceCenter[1] / count;
            dualGrid.vertices[i].coords[2] = faceCenter[2] / count;
        }
        // set edges
        for (var i = this.edges.length - 1; i >= 0; i--) {
            var edge = this.edges[i];
            dualGrid.edges[i].faceA = edge.firstVertexA;
            dualGrid.edges[i].faceB = edge.firstVertexB;
            // faces are swapped from verts
            dualGrid.edges[i].firstVertexA = edge.faceB;
            dualGrid.edges[i].firstVertexB = edge.faceA;
        }
        // set prev and next
        for (var faceIndex = dualGrid.faces.length - 1; faceIndex >= 0; faceIndex--) {
            var face = dualGrid.faces[faceIndex];
            for (var faceEdgeIndex = face.edges.length - 1; faceEdgeIndex >= 0; faceEdgeIndex--) {
                var edgeIndex = face.edges[faceEdgeIndex];
                var theEdge = dualGrid.edges[edgeIndex];
                if (theEdge.faceA == faceIndex) {
                    if (faceEdgeIndex == 0) {
                        dualGrid.edges[edgeIndex].prevA = face.edges[face.edges.length - 1];
                        dualGrid.edges[edgeIndex].nextA = face.edges[1];
                    }
                    else if (faceEdgeIndex == face.edges.length - 1) {
                        dualGrid.edges[edgeIndex].prevA = face.edges[face.edges.length - 1];
                        dualGrid.edges[edgeIndex].nextA = face.edges[0];
                    }
                    else {
                        dualGrid.edges[edgeIndex].prevA = face.edges[faceEdgeIndex - 1];
                        dualGrid.edges[edgeIndex].nextA = face.edges[faceEdgeIndex + 1];
                    }
                }
                if (theEdge.faceB == faceIndex) {
                    if (faceEdgeIndex == 0) {
                        dualGrid.edges[edgeIndex].prevB = face.edges[face.edges.length - 1];
                        dualGrid.edges[edgeIndex].nextB = face.edges[1];
                    }
                    else if (faceEdgeIndex == face.edges.length - 1) {
                        dualGrid.edges[edgeIndex].prevB = face.edges[face.edges.length - 1];
                        dualGrid.edges[edgeIndex].nextB = face.edges[0];
                    }
                    else {
                        dualGrid.edges[edgeIndex].prevB = face.edges[faceEdgeIndex - 1];
                        dualGrid.edges[edgeIndex].nextB = face.edges[faceEdgeIndex + 1];
                    }
                }
            }
        }
        return [dualGrid, null];
    };
    /************* Subdivision helpers ****************/
    // rename vertexInexAtClockwiseIndexOnFaceToSubdivide
    WingedGrid.prototype.vertexIndexAtClockwiseIndexOnOldFace = function (faceIndex, edgeInFaceIndex, clockwiseVertexIndex, edgeSubdivisions) {
        var edgeIndex = this.faces[faceIndex].edges[edgeInFaceIndex];
        var edge = this.edges[edgeIndex];
        if (edge.faceA == faceIndex) {
            return this.vertices.length + edgeIndex * edgeSubdivisions + clockwiseVertexIndex;
        }
        if (edge.faceB == faceIndex) {
            return this.vertices.length + edgeIndex * edgeSubdivisions + edgeSubdivisions - 1 - clockwiseVertexIndex;
        }
        return -1;
    };
    WingedGrid.prototype.edgeIndexAtClockwiseIndexOnOldFace = function (faceIndex, edgeInFaceIndex, clockwiseEdgeIndex, edgeSubdivisions) {
        var oldEdgeIndex = this.faces[faceIndex].edges[edgeInFaceIndex];
        var oldEdge = this.edges[oldEdgeIndex];
        if (oldEdge.faceA == faceIndex) {
            return oldEdgeIndex * (edgeSubdivisions + 1) + clockwiseEdgeIndex;
        }
        if (oldEdge.faceB == faceIndex) {
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
        var faceCount = this.faces.length * ((edgeSubdivisions + 2) * (edgeSubdivisions + 1) / 2 + (edgeSubdivisions + 1) * edgeSubdivisions / 2);
        var dividedGrid = new WingedGrid(faceCount, 
        // since each face 'owns' 1/2 of three edges, there are 1.5 times as
        //  many edges as faces
        3 * faceCount / 2, 
        // Euler-somebody or other gives us the vertex count of
        //  1/2*(face count) + 2
        faceCount / 2 + 2);
        // Invalidate all values that will be set later, useful for error checking
        for (var i_1 = 0; i_1 < faceCount; i_1++) {
            dividedGrid.faces[i_1].edges = [-1, -1, -1];
        }
        for (var i_2 = 0; i_2 < 3 * faceCount / 2; i_2++) {
            dividedGrid.edges[i_2].faceA = -1;
            dividedGrid.edges[i_2].faceB = -1;
            dividedGrid.edges[i_2].firstVertexA = -1;
            dividedGrid.edges[i_2].firstVertexB = -1;
            dividedGrid.edges[i_2].nextA = -1;
            dividedGrid.edges[i_2].nextB = -1;
            dividedGrid.edges[i_2].prevA = -1;
            dividedGrid.edges[i_2].prevB = -1;
        }
        // vertices corrisponding with old ones will have the same number of
        //  associated edges to preserve the Euler characteristic ( =2 for S2)
        for (var i_3 = 0; i_3 < this.vertices.length; i_3++) {
            dividedGrid.vertices[i_3].edges = new Array(this.vertices[i_3].edges.length);
            for (var j = 0; j < this.vertices[i_3].edges.length; j++) {
                dividedGrid.vertices[i_3].edges[j] = -1;
            }
            // set the coords
            dividedGrid.vertices[i_3].coords = [Infinity, Infinity, Infinity];
        }
        // the way we divide the faces creates six accociated edges for the remaining
        //  vertecies
        for (var i_4 = this.vertices.length; i_4 < faceCount / 2 + 2; i_4++) {
            dividedGrid.vertices[i_4].edges = [-1, -1, -1, -1, -1, -1];
            // set the coords
            dividedGrid.vertices[i_4].coords = [Infinity, Infinity, Infinity];
        }
        /***************** Subdivide the grid ****************/
        /******************* EDGE SUBDIVISION ***********************/
        /****** Old Edge Subdivision goes in the first section of the new array, ordered by edge *******/
        var origVertexCount = this.vertices.length;
        var origEdgeCount = this.edges.length;
        for (var i_5 = this.edges.length - 1; i_5 >= 0; i_5--) {
            var edge_1 = this.edges[i_5];
            // first edge has origional vertex
            dividedGrid.edges[i_5 * (edgeSubdivisions + 1)].firstVertexA = edge_1.firstVertexA;
            dividedGrid.edges[i_5 * (edgeSubdivisions + 1)].firstVertexB = origVertexCount + i_5 * edgeSubdivisions;
            for (var j = 1; j < edgeSubdivisions; j++) {
                // set the edge vertex indecies
                dividedGrid.edges[i_5 * (edgeSubdivisions + 1) + j].firstVertexA = origVertexCount + i_5 * edgeSubdivisions + j - 1;
                dividedGrid.edges[i_5 * (edgeSubdivisions + 1) + j].firstVertexB = origVertexCount + i_5 * edgeSubdivisions + j;
            }
            // connect last new edge
            dividedGrid.edges[i_5 * (edgeSubdivisions + 1) + edgeSubdivisions].firstVertexA = origVertexCount + i_5 * edgeSubdivisions + edgeSubdivisions - 1;
            dividedGrid.edges[i_5 * (edgeSubdivisions + 1) + edgeSubdivisions].firstVertexB = edge_1.firstVertexB;
        }
        /********* Edges created interior to old faces go in the second section, ordered by face. ****/
        for (var faceIndex_1 = this.faces.length - 1; faceIndex_1 >= 0; faceIndex_1--) {
            var oldFace = this.faces[faceIndex_1];
            // vertex offset for vertices interior to the face
            var vertexOffset_1 = origVertexCount + origEdgeCount * edgeSubdivisions + (edgeSubdivisions * (edgeSubdivisions - 1) / 2) * faceIndex_1;
            var edgeOffset = (edgeSubdivisions + 1) * origEdgeCount + 3 * (edgeSubdivisions) * (edgeSubdivisions + 1) / 2 * faceIndex_1;
            if (edgeSubdivisions == 1) {
                dividedGrid.edges[edgeOffset + 0].firstVertexA = origVertexCount + oldFace.edges[0];
                dividedGrid.edges[edgeOffset + 0].firstVertexB = origVertexCount + oldFace.edges[1];
                dividedGrid.edges[edgeOffset + 1].firstVertexA = origVertexCount + oldFace.edges[0];
                dividedGrid.edges[edgeOffset + 1].firstVertexB = origVertexCount + oldFace.edges[2];
                dividedGrid.edges[edgeOffset + 2].firstVertexA = origVertexCount + oldFace.edges[1];
                dividedGrid.edges[edgeOffset + 2].firstVertexB = origVertexCount + oldFace.edges[2];
            }
            else {
                // first row
                dividedGrid.edges[edgeOffset + 0].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + 0].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, 0, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + 1].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + 1].firstVertexB = vertexOffset_1 + 0;
                dividedGrid.edges[edgeOffset + 2].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, 0, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + 2].firstVertexB = vertexOffset_1 + 0;
                // middle rows
                var rowOffset;
                var i = 0;
                for (i = 1; i < edgeSubdivisions - 1; i++) {
                    rowOffset = i * (i + 1) * 3 / 2;
                    // first border
                    dividedGrid.edges[edgeOffset + rowOffset + 0].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1 - i, edgeSubdivisions);
                    dividedGrid.edges[edgeOffset + rowOffset + 0].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2);
                    dividedGrid.edges[edgeOffset + rowOffset + 1].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, edgeSubdivisions - 1 - i, edgeSubdivisions);
                    dividedGrid.edges[edgeOffset + rowOffset + 1].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + i;
                    dividedGrid.edges[edgeOffset + rowOffset + 2].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2);
                    dividedGrid.edges[edgeOffset + rowOffset + 2].firstVertexB = vertexOffset_1 + (i * (i - 1)) / 2 + i;
                    // interior of face
                    for (var j = 1; j < i; j++) {
                        dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 0].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                        dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 0].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + j;
                        dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 1].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                        dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 1].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + i + j;
                        dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 2].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + j;
                        dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 2].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + i + j;
                    }
                    // second border
                    dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 0].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                    dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 0].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                    dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 1].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                    dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 1].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + 2 * i;
                    dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 2].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                    dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 2].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + 2 * i;
                }
                // last row
                i = edgeSubdivisions - 1; // should already be set, but just incase
                rowOffset = i * (i + 1) * 3 / 2;
                // border 1
                dividedGrid.edges[edgeOffset + rowOffset + 0].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, 0, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + rowOffset + 0].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2);
                dividedGrid.edges[edgeOffset + rowOffset + 1].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 0, 0, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + rowOffset + 1].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - 1, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + rowOffset + 2].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2);
                dividedGrid.edges[edgeOffset + rowOffset + 2].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - 1, edgeSubdivisions);
                // middle
                for (var j = 1; j < i; j++) {
                    dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 0].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                    dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 0].firstVertexB = vertexOffset_1 + (i * (i - 1) / 2) + j;
                    dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 1].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + j - 1;
                    dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 1].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - j - 1, edgeSubdivisions);
                    dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 2].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + j;
                    dividedGrid.edges[edgeOffset + rowOffset + j * 3 + 2].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, edgeSubdivisions - j - 1, edgeSubdivisions);
                }
                // border 2
                dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 0].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 0].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 1].firstVertexA = vertexOffset_1 + (i * (i - 1) / 2) + i - 1;
                dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 1].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, 0, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 2].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 1, i, edgeSubdivisions);
                dividedGrid.edges[edgeOffset + rowOffset + i * 3 + 2].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_1, 2, 0, edgeSubdivisions);
            }
        }
        // create the faces, update the edges
        /******************* FACE SUBDIVISION ***********************/
        for (var faceIndex_2 = this.faces.length - 1; faceIndex_2 >= 0; faceIndex_2--) {
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
            dividedGrid.faces[indexStart + 0].edges =
                [this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 0, edgeSubdivisions, edgeSubdivisions),
                    this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 1, 0, edgeSubdivisions),
                    edgeOffset];
            // loop through the middle section of faces
            // edges grow by 1/2*i*(i-1)*3
            for (var i_6 = 1; i_6 < edgeSubdivisions; i_6++) {
                var rowIndexStart = i_6 * (i_6 + 1) / 2 + i_6 * (i_6 - 1) / 2;
                // edge
                dividedGrid.faces[indexStart + rowIndexStart + 0].edges =
                    [this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 0, edgeSubdivisions - i_6, edgeSubdivisions),
                        edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + 1,
                        edgeOffset + i_6 * (i_6 + 1) * 3 / 2];
                // middle
                // up to one less than next index start
                for (var j = 1; j < (i_6 + 1) * (i_6 + 2) / 2 + (i_6 + 1) * i_6 / 2 - rowIndexStart - 1; j++) {
                    if (j % 2 == 1) {
                        dividedGrid.faces[indexStart + rowIndexStart + j].edges =
                            [edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + (j - 1) * 3 / 2,
                                edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + (j - 1) * 3 / 2 + 2,
                                edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + (j - 1) * 3 / 2 + 1];
                    }
                    else {
                        dividedGrid.faces[indexStart + rowIndexStart + j].edges =
                            [edgeOffset + i_6 * (i_6 + 1) * 3 / 2 + j * 3 / 2,
                                edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + (j - 2) * 3 / 2 + 2,
                                edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + j * 3 / 2 + 1];
                    }
                }
                // edge
                dividedGrid.faces[indexStart + (i_6 + 1) * (i_6 + 2) / 2 + (i_6 + 1) * i_6 / 2 - 1].edges =
                    [edgeOffset + i_6 * (i_6 + 1) * 3 / 2 + ((i_6 + 1) * (i_6 + 2) / 2 + (i_6 + 1) * i_6 / 2 - rowIndexStart - 1) * 3 / 2,
                        edgeOffset + i_6 * (i_6 - 1) * 3 / 2 + ((i_6 + 1) * (i_6 + 2) / 2 + (i_6 + 1) * i_6 / 2 - rowIndexStart - 3) * 3 / 2 + 2,
                        this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 1, i_6, edgeSubdivisions)];
            }
            var rowIndexStart = edgeSubdivisions * (edgeSubdivisions + 1) / 2 + edgeSubdivisions * (edgeSubdivisions - 1) / 2;
            // bottom corner 1
            dividedGrid.faces[indexStart + rowIndexStart + 0].edges =
                [this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 0, 0, edgeSubdivisions),
                    edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + 1,
                    this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 2, edgeSubdivisions, edgeSubdivisions)];
            // bottom edge
            for (var j = 1; j < subFaceCount - rowIndexStart - 1; j++) {
                if (j % 2 == 1) {
                    dividedGrid.faces[indexStart + rowIndexStart + j].edges =
                        [edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2,
                            edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2 + 2,
                            edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2 + 1];
                }
                else {
                    dividedGrid.faces[indexStart + rowIndexStart + j].edges =
                        [this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 2, edgeSubdivisions - j / 2, edgeSubdivisions),
                            edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 2) * 3 / 2 + 2,
                            edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + j * 3 / 2 + 1];
                }
            }
            // bottom corner 2
            dividedGrid.faces[indexStart + subFaceCount - 1].edges =
                [this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 2, 0, edgeSubdivisions),
                    edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (subFaceCount - rowIndexStart - 3) * 3 / 2 + 2,
                    this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex_2, 1, edgeSubdivisions, edgeSubdivisions)];
        }
        // set edge faces from the previously build edge arrays
        for (var faceIndex = dividedGrid.faces.length - 1; faceIndex >= 0; faceIndex--) {
            var edges = dividedGrid.faces[faceIndex].edges;
            var edgesLength = edges.length;
            var thisEdge, nextEdge, prevEdge;
            /*for (let i = 0; i < edges.length-1; i++) {
                const prevIndex: number = (i + edges.length - 1) % edges.length;
                const nextIndex: number = (i + 1) % edges.length;


                prevEdge = dividedGrid.edges[edges[prevIndex]];
                thisEdge = dividedGrid.edges[edges[i]];
                nextEdge = dividedGrid.edges[edges[nextIndex]];
                // test vertices
                if (thisEdge.firstVertexA == prevEdge.firstVertexA || thisEdge.firstVertexA == prevEdge.firstVertexB) {
                    // check the next edge also matches
                    if (thisEdge.firstVertexB == nextEdge.firstVertexA || thisEdge.firstVertexB == nextEdge.firstVertexB) {
                        if (thisEdge.faceA == -1) {
                            dividedGrid.edges[edges[i]].faceA = faceIndex;
                            dividedGrid.edges[edges[i]].prevA = edges[prevIndex];
                            dividedGrid.edges[edges[i]].nextA = edges[nextIndex];
                        } else {
                            console.log("For face "+faceIndex+". Face A has already been set for edge: "+thisEdge+" With edge set: "+edges);
                        }

                    } else {
                        console.log("For face "+faceIndex+". Previous edge matches, but next edge doesn't share correct vertex!");
                    }
                } else if (thisEdge.firstVertexB == prevEdge.firstVertexA || thisEdge.firstVertexB == prevEdge.firstVertexB) {
                    // check the next edge also matches
                    if (thisEdge.firstVertexA == nextEdge.firstVertexA || thisEdge.firstVertexA == nextEdge.firstVertexB) {
                        if (thisEdge.faceB == -1) {
                            dividedGrid.edges[edges[i]].faceB = faceIndex;
                            dividedGrid.edges[edges[i]].prevB = edges[prevIndex];
                            dividedGrid.edges[edges[i]].nextB = edges[nextIndex];
                        } else {
                            console.log("For face "+faceIndex+". Face B has already been set for edge: "+thisEdge+" With edge set: "+edges);
                        }

                    } else {
                        console.log("For face "+faceIndex+". Previous edge matches, but next edge doesn't share correct vertex!");
                    }
                } else {
                    console.log("For face "+faceIndex+". Edges Don't share a vertex!");
                }
            }//*/
            // ****** Old method, keep until above is known to work
            var prevEdge_1 = dividedGrid.edges[edges[edgesLength - 1]];
            var thisEdge_1 = dividedGrid.edges[edges[0]];
            var nextEdge_1 = dividedGrid.edges[edges[1]];
            // test vertices
            if (thisEdge_1.firstVertexA == prevEdge_1.firstVertexA || thisEdge_1.firstVertexA == prevEdge_1.firstVertexB) {
                // check the next edge also matches and face A is not set
                if (thisEdge_1.firstVertexB == nextEdge_1.firstVertexA || thisEdge_1.firstVertexB == nextEdge_1.firstVertexB) {
                    if (thisEdge_1.faceA == -1) {
                        dividedGrid.edges[edges[0]].faceA = faceIndex;
                        dividedGrid.edges[edges[0]].prevA = edges[edgesLength - 1];
                        dividedGrid.edges[edges[0]].nextA = edges[1];
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge_1 + " With edge set: " + edges);
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                }
            }
            else if (thisEdge_1.firstVertexB == prevEdge_1.firstVertexA || thisEdge_1.firstVertexB == prevEdge_1.firstVertexB) {
                // check the next edge also matches and face B is not set
                if (thisEdge_1.firstVertexA == nextEdge_1.firstVertexA || thisEdge_1.firstVertexA == nextEdge_1.firstVertexB) {
                    if (thisEdge_1.faceB == -1) {
                        dividedGrid.edges[edges[0]].faceB = faceIndex;
                        dividedGrid.edges[edges[0]].prevB = edges[edgesLength - 1];
                        dividedGrid.edges[edges[0]].nextB = edges[1];
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge_1 + " With edge set: " + edges);
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
            for (var i_7 = 1; i_7 < edgesLength - 1; i_7++) {
                prevEdge_1 = dividedGrid.edges[edges[i_7 - 1]];
                thisEdge_1 = dividedGrid.edges[edges[i_7]];
                nextEdge_1 = dividedGrid.edges[edges[i_7 + 1]];
                // test vertcies
                if (thisEdge_1.firstVertexA == prevEdge_1.firstVertexA || thisEdge_1.firstVertexA == prevEdge_1.firstVertexB) {
                    // check the next edge also matches
                    if (thisEdge_1.firstVertexB == nextEdge_1.firstVertexA || thisEdge_1.firstVertexB == nextEdge_1.firstVertexB) {
                        if (thisEdge_1.faceA == -1) {
                            dividedGrid.edges[edges[i_7]].faceA = faceIndex;
                            dividedGrid.edges[edges[i_7]].prevA = edges[i_7 - 1];
                            dividedGrid.edges[edges[i_7]].nextA = edges[i_7 + 1];
                        }
                        else {
                            console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge_1 + " With edge set: " + edges);
                        }
                    }
                    else {
                        console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                    }
                }
                else if (thisEdge_1.firstVertexB == prevEdge_1.firstVertexA || thisEdge_1.firstVertexB == prevEdge_1.firstVertexB) {
                    // check the next edge also matches
                    if (thisEdge_1.firstVertexA == nextEdge_1.firstVertexA || thisEdge_1.firstVertexA == nextEdge_1.firstVertexB) {
                        if (thisEdge_1.faceB == -1) {
                            dividedGrid.edges[edges[i_7]].faceB = faceIndex;
                            dividedGrid.edges[edges[i_7]].prevB = edges[i_7 - 1];
                            dividedGrid.edges[edges[i_7]].nextB = edges[i_7 + 1];
                        }
                        else {
                            console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge_1 + " With edge set: " + edges);
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
            prevEdge_1 = dividedGrid.edges[edges[edgesLength - 2]];
            thisEdge_1 = dividedGrid.edges[edges[edgesLength - 1]];
            nextEdge_1 = dividedGrid.edges[edges[0]];
            // test vertices
            if (thisEdge_1.firstVertexA == prevEdge_1.firstVertexA || thisEdge_1.firstVertexA == prevEdge_1.firstVertexB) {
                // check the next edge also matches
                if (thisEdge_1.firstVertexB == nextEdge_1.firstVertexA || thisEdge_1.firstVertexB == nextEdge_1.firstVertexB) {
                    if (thisEdge_1.faceA == -1) {
                        dividedGrid.edges[edges[edgesLength - 1]].faceA = faceIndex;
                        dividedGrid.edges[edges[edgesLength - 1]].prevA = edges[edgesLength - 2];
                        dividedGrid.edges[edges[edgesLength - 1]].nextA = edges[0];
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge_1 + " With edge set: " + edges);
                    }
                }
                else {
                    console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
                }
            }
            else if (thisEdge_1.firstVertexB == prevEdge_1.firstVertexA || thisEdge_1.firstVertexB == prevEdge_1.firstVertexB) {
                // check the next edge also matches
                if (thisEdge_1.firstVertexA == nextEdge_1.firstVertexA || thisEdge_1.firstVertexA == nextEdge_1.firstVertexB) {
                    if (thisEdge_1.faceB == -1) {
                        dividedGrid.edges[edges[edgesLength - 1]].faceB = faceIndex;
                        dividedGrid.edges[edges[edgesLength - 1]].prevB = edges[edgesLength - 2];
                        dividedGrid.edges[edges[edgesLength - 1]].nextB = edges[0];
                    }
                    else {
                        console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge_1 + " With edge set: " + edges);
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
        for (var i_8 = this.vertices.length - 1; i_8 >= 0; i_8--) {
            var vertex = this.vertices[i_8];
            dividedGrid.vertices[i_8].coords[0] = vertex.coords[0];
            dividedGrid.vertices[i_8].coords[1] = vertex.coords[1];
            dividedGrid.vertices[i_8].coords[2] = vertex.coords[2];
        }
        // subdivide along each edge
        for (var i_9 = this.edges.length - 1; i_9 >= 0; i_9--) {
            var edge_2 = this.edges[i_9];
            var firstVertex_1 = this.vertices[edge_2.firstVertexA];
            var secondVertex_1 = this.vertices[edge_2.firstVertexB];
            // angle to subdivide
            var angleToSubdivide_1 = void 0;
            angleToSubdivide_1 = vectorAngle(firstVertex_1.coords, secondVertex_1.coords);
            // angle between origin, first vertex, and second vertex
            var vectorA_1 = [-1, -1, -1], vectorB_1 = [-1, -1, -1];
            vectorA_1[0] = -1 * firstVertex_1.coords[0];
            vectorA_1[1] = -1 * firstVertex_1.coords[1];
            vectorA_1[2] = -1 * firstVertex_1.coords[2];
            vectorB_1[0] = secondVertex_1.coords[0] - firstVertex_1.coords[0];
            vectorB_1[1] = secondVertex_1.coords[1] - firstVertex_1.coords[1];
            vectorB_1[2] = secondVertex_1.coords[2] - firstVertex_1.coords[2];
            var cornerAngle_1 = void 0;
            cornerAngle_1 = vectorAngle(vectorA_1, vectorB_1);
            // unit vector from first to second vertex
            var stepDirection_1 = [-1, -1, -1];
            stepDirection_1[0] = vectorB_1[0] / vectorLength(vectorB_1);
            stepDirection_1[1] = vectorB_1[1] / vectorLength(vectorB_1);
            stepDirection_1[2] = vectorB_1[2] / vectorLength(vectorB_1);
            // origional radius of the
            var sphereRadius = vectorLength(firstVertex_1.coords);
            var divisionLength_1 = void 0;
            for (var j = 0; j < edgeSubdivisions; j++) {
                // find the new vertex position and create the vertex
                // but don't correct it's length yet
                divisionLength_1 = Math.sin(angleToSubdivide_1 * ((j + 1) / (edgeSubdivisions + 1))) * sphereRadius / Math.sin(Math.PI - cornerAngle_1 - angleToSubdivide_1 * ((j + 1) / (edgeSubdivisions + 1)));
                dividedGrid.vertices[origVertexCount + i_9 * edgeSubdivisions + j].coords[0] = firstVertex_1.coords[0] + stepDirection_1[0] * divisionLength_1;
                dividedGrid.vertices[origVertexCount + i_9 * edgeSubdivisions + j].coords[1] = firstVertex_1.coords[1] + stepDirection_1[1] * divisionLength_1;
                dividedGrid.vertices[origVertexCount + i_9 * edgeSubdivisions + j].coords[2] = firstVertex_1.coords[2] + stepDirection_1[2] * divisionLength_1;
            }
        }
        // subdivide face interior
        for (var faceIndex_3 = 0; faceIndex_3 < this.faces.length; faceIndex_3++) {
            var vertexOffset = origVertexCount + origEdgeCount * edgeSubdivisions + edgeSubdivisions * (edgeSubdivisions - 1) / 2 * faceIndex_3;
            // only if we have more than one division
            if (edgeSubdivisions > 1) {
                for (var i_10 = 0; i_10 < edgeSubdivisions - 1; i_10++) {
                    var firstVertex = dividedGrid.vertices[this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_3, 0, edgeSubdivisions - 2 - i_10, edgeSubdivisions)];
                    var secondVertex = dividedGrid.vertices[this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex_3, 1, 1 + i_10, edgeSubdivisions)];
                    // angle to subdivide
                    var angleToSubdivide;
                    angleToSubdivide = vectorAngle(firstVertex.coords, secondVertex.coords);
                    // angle between origin, first vertex, and second vertex
                    var vectorA = [-1, -1, -1], vectorB = [-1, -1, -1];
                    vectorA[0] = -1 * firstVertex.coords[0];
                    vectorA[1] = -1 * firstVertex.coords[1];
                    vectorA[2] = -1 * firstVertex.coords[2];
                    vectorB[0] = secondVertex.coords[0] - firstVertex.coords[0];
                    vectorB[1] = secondVertex.coords[1] - firstVertex.coords[1];
                    vectorB[2] = secondVertex.coords[2] - firstVertex.coords[2];
                    var cornerAngle;
                    cornerAngle = vectorAngle(vectorA, vectorB);
                    // unit vector from first to second vertex
                    var stepDirection = [-1, -1, -1];
                    stepDirection[0] = vectorB[0] / vectorLength(vectorB);
                    stepDirection[1] = vectorB[1] / vectorLength(vectorB);
                    stepDirection[2] = vectorB[2] / vectorLength(vectorB);
                    var firstVectorLength = vectorLength(firstVertex.coords);
                    var divisionLength = 0;
                    for (var j = 0; j < i_10 + 1; j++) {
                        divisionLength = Math.sin(angleToSubdivide * ((j + 1) / (i_10 + 2))) * firstVectorLength / Math.sin(Math.PI - cornerAngle - angleToSubdivide * ((j + 1) / (i_10 + 2)));
                        if (vertexOffset + (i_10 * (i_10 + 1) / 2) + j >= dividedGrid.vertices.length) {
                            console.log("breakpoint");
                        }
                        dividedGrid.vertices[vertexOffset + (i_10 * (i_10 + 1) / 2) + j].coords[0] = firstVertex.coords[0] + stepDirection[0] * divisionLength;
                        dividedGrid.vertices[vertexOffset + (i_10 * (i_10 + 1) / 2) + j].coords[1] = firstVertex.coords[1] + stepDirection[1] * divisionLength;
                        dividedGrid.vertices[vertexOffset + (i_10 * (i_10 + 1) / 2) + j].coords[2] = firstVertex.coords[2] + stepDirection[2] * divisionLength;
                    }
                }
            }
        }
        // set vertex edge array
        // loop through edges so we only have to touch each one once
        for (var index = dividedGrid.edges.length - 1; index >= 0; index--) {
            var edge = dividedGrid.edges[index];
            var nextEdgeIndex = -1;
            var nextEdge;
            var theVertexIndex = edge.firstVertexA;
            var theVertex = dividedGrid.vertices[theVertexIndex];
            var error;
            if (theVertex.edges[0] == -1) {
                _a = edge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _a[0], error = _a[1];
                nextEdge = dividedGrid.edges[nextEdgeIndex];
                theVertex.edges[0] = index;
                var i = 1;
                while (index != nextEdgeIndex) {
                    theVertex.edges[i] = nextEdgeIndex;
                    _b = nextEdge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _b[0], error = _b[1];
                    nextEdge = dividedGrid.edges[nextEdgeIndex];
                    i = i + 1;
                }
            }
            theVertexIndex = edge.firstVertexB;
            theVertex = dividedGrid.vertices[theVertexIndex];
            if (theVertex.edges[0] == -1) {
                _c = edge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _c[0], error = _c[1];
                nextEdge = dividedGrid.edges[nextEdgeIndex];
                theVertex.edges[0] = index;
                var i = 1;
                while (index != nextEdgeIndex) {
                    theVertex.edges[i] = nextEdgeIndex;
                    _d = nextEdge.nextEdgeForVertex(theVertexIndex), nextEdgeIndex = _d[0], error = _d[1];
                    nextEdge = dividedGrid.edges[nextEdgeIndex];
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
    var icosahedron = new WingedGrid(0, 0, 0); // arrays will be set manually
    icosahedron.faces = [
        // cap 1
        { edges: [9, 4, 0] } // 0
        ,
        { edges: [1, 5, 0] } // 1
        ,
        { edges: [2, 6, 1] } // 2
        ,
        { edges: [3, 7, 2] } // 3
        ,
        { edges: [4, 8, 3] } // 4
        ,
        { edges: [14, 19, 10] } // 5
        ,
        { edges: [15, 11, 10] } // 6
        ,
        { edges: [16, 12, 11] } // 7
        ,
        { edges: [17, 13, 12] } // 8
        ,
        { edges: [18, 14, 13] } // 9
        ,
        { edges: [29, 28, 9] } // 10
        ,
        { edges: [29, 20, 17] } // 11
        ,
        { edges: [21, 20, 5] } // 12
        ,
        { edges: [21, 22, 18] } // 13
        ,
        { edges: [23, 22, 6] } // 14
        ,
        { edges: [23, 24, 19] } // 15
        ,
        { edges: [25, 24, 7] } // 16
        ,
        { edges: [25, 26, 15] } // 17
        ,
        { edges: [27, 26, 8] } // 18
        ,
        { edges: [27, 28, 16] } // 19
    ];
    var goldenRatio = 1.61803398875;
    icosahedron.vertices = [
        // Edges are duplicate info, could be done pragmatically after constucting
        // edges
        // Must be in counter-clockwise order
        // y-z plane rectangle
        { coords: [0, 1, goldenRatio],
            edges: [4, 3, 2, 1, 0],
            vertexNeighbors: [] } // 0
        ,
        { coords: [0, 1, -goldenRatio],
            edges: [19, 24, 25, 15, 10],
            vertexNeighbors: [] } // 1
        ,
        { coords: [0, -1, goldenRatio],
            edges: [5, 20, 29, 9, 0],
            vertexNeighbors: [] } // 2
        ,
        { coords: [0, -1, -goldenRatio],
            edges: [11, 12, 13, 14, 10],
            vertexNeighbors: [] } // 3
        ,
        { coords: [-goldenRatio, 0, 1],
            edges: [6, 22, 21, 5, 1],
            vertexNeighbors: [] } // 4
        ,
        { coords: [goldenRatio, 0, 1],
            edges: [9, 28, 27, 8, 4],
            vertexNeighbors: [] } // 5
        ,
        { coords: [-goldenRatio, 0, -1],
            edges: [18, 22, 23, 19, 14],
            vertexNeighbors: [] } // 6
        ,
        { coords: [goldenRatio, 0, -1],
            edges: [15, 26, 27, 16, 11],
            vertexNeighbors: [] } // 7
        ,
        { coords: [-1, goldenRatio, 0],
            edges: [7, 24, 23, 6, 2],
            vertexNeighbors: [] } // 8
        ,
        { coords: [-1, -goldenRatio, 0],
            edges: [17, 20, 21, 18, 13],
            vertexNeighbors: [] } // 9
        ,
        { coords: [1, goldenRatio, 0],
            edges: [8, 26, 25, 7, 3],
            vertexNeighbors: [] } // 10
        ,
        { coords: [1, -goldenRatio, 0],
            edges: [16, 28, 29, 17, 12],
            vertexNeighbors: [] } // 11
    ];
    icosahedron.edges = [
        // cap 1 around vertex 0
        // 5 spokes starting between face 0 and 1,
        // clockwise around starting with short edge
        // of the rectangle (y-z)
        {
            firstVertexA: 0, firstVertexB: 2,
            faceA: 0, faceB: 1,
            prevA: 4, nextA: 9,
            prevB: 5, nextB: 1
        },
        {
            firstVertexA: 0, firstVertexB: 4,
            faceA: 1, faceB: 2,
            prevA: 0, nextA: 5,
            prevB: 6, nextB: 2
        },
        {
            firstVertexA: 0, firstVertexB: 8,
            faceA: 2, faceB: 3,
            prevA: 1, nextA: 6,
            prevB: 7, nextB: 3
        },
        {
            firstVertexA: 0, firstVertexB: 10,
            faceA: 3, faceB: 4,
            prevA: 2, nextA: 7,
            prevB: 8, nextB: 4
        },
        {
            firstVertexA: 0, firstVertexB: 5,
            faceA: 4, faceB: 0,
            prevA: 3, nextA: 8,
            prevB: 9, nextB: 0
        },
        // ring of 5 around the base of the cap
        {
            firstVertexA: 4, firstVertexB: 2,
            faceA: 1, faceB: 12,
            prevA: 1, nextA: 0,
            prevB: 20, nextB: 21
        },
        {
            firstVertexA: 8, firstVertexB: 4,
            faceA: 2, faceB: 14,
            prevA: 2, nextA: 1,
            prevB: 22, nextB: 23
        },
        {
            firstVertexA: 10, firstVertexB: 8,
            faceA: 3, faceB: 16,
            prevA: 3, nextA: 2,
            prevB: 24, nextB: 25
        },
        {
            firstVertexA: 5, firstVertexB: 10,
            faceA: 4, faceB: 18,
            prevA: 4, nextA: 3,
            prevB: 26, nextB: 27
        },
        {
            firstVertexA: 2, firstVertexB: 5,
            faceA: 0, faceB: 10,
            prevA: 0, nextA: 4,
            prevB: 28, nextB: 29
        },
        // cap 2 around vertex 3
        // 5 spokes starting between face 5 and 6,
        // counter-clockwise from short edge (y-z) rectangle
        {
            firstVertexA: 1, firstVertexB: 3,
            faceA: 5, faceB: 6,
            prevA: 19, nextA: 14,
            prevB: 11, nextB: 15
        },
        {
            firstVertexA: 7, firstVertexB: 3,
            faceA: 6, faceB: 7,
            prevA: 15, nextA: 10,
            prevB: 12, nextB: 16
        },
        {
            firstVertexA: 11, firstVertexB: 3,
            faceA: 7, faceB: 8,
            prevA: 16, nextA: 11,
            prevB: 13, nextB: 17
        },
        {
            firstVertexA: 9, firstVertexB: 3,
            faceA: 8, faceB: 9,
            prevA: 17, nextA: 12,
            prevB: 14, nextB: 18
        },
        {
            firstVertexA: 6, firstVertexB: 3,
            faceA: 9, faceB: 5,
            prevA: 18, nextA: 13,
            prevB: 10, nextB: 19
        },
        // ring of 5 around the base of cap 2
        {
            firstVertexA: 1, firstVertexB: 7,
            faceA: 6, faceB: 17,
            prevA: 10, nextA: 11,
            prevB: 26, nextB: 25
        },
        {
            firstVertexA: 7, firstVertexB: 11,
            faceA: 7, faceB: 19,
            prevA: 11, nextA: 12,
            prevB: 28, nextB: 27
        },
        {
            firstVertexA: 11, firstVertexB: 9,
            faceA: 8, faceB: 11,
            prevA: 12, nextA: 13,
            prevB: 20, nextB: 29
        },
        {
            firstVertexA: 9, firstVertexB: 6,
            faceA: 9, faceB: 13,
            prevA: 13, nextA: 14,
            prevB: 22, nextB: 21
        },
        {
            firstVertexA: 6, firstVertexB: 1,
            faceA: 5, faceB: 15,
            prevA: 14, nextA: 10,
            prevB: 24, nextB: 23
        },
        // zig-zag down the middle
        // 10 triangles, 10 new edges
        // starting clockwise from end of edge 0
        {
            firstVertexA: 2, firstVertexB: 9,
            faceA: 11, faceB: 12,
            prevA: 29, nextA: 17,
            prevB: 21, nextB: 5
        },
        {
            firstVertexA: 4, firstVertexB: 9,
            faceA: 12, faceB: 13,
            prevA: 5, nextA: 20,
            prevB: 18, nextB: 22
        },
        {
            firstVertexA: 4, firstVertexB: 6,
            faceA: 13, faceB: 14,
            prevA: 21, nextA: 18,
            prevB: 23, nextB: 6
        },
        {
            firstVertexA: 8, firstVertexB: 6,
            faceA: 14, faceB: 15,
            prevA: 6, nextA: 22,
            prevB: 19, nextB: 24
        },
        {
            firstVertexA: 8, firstVertexB: 1,
            faceA: 15, faceB: 16,
            prevA: 23, nextA: 19,
            prevB: 25, nextB: 7
        },
        {
            firstVertexA: 10, firstVertexB: 1,
            faceA: 16, faceB: 17,
            prevA: 7, nextA: 24,
            prevB: 15, nextB: 26
        },
        {
            firstVertexA: 10, firstVertexB: 7,
            faceA: 17, faceB: 18,
            prevA: 25, nextA: 15,
            prevB: 27, nextB: 8
        },
        {
            firstVertexA: 5, firstVertexB: 7,
            faceA: 18, faceB: 19,
            prevA: 8, nextA: 26,
            prevB: 16, nextB: 28
        },
        {
            firstVertexA: 5, firstVertexB: 11,
            faceA: 19, faceB: 10,
            prevA: 27, nextA: 16,
            prevB: 29, nextB: 9
        },
        {
            firstVertexA: 2, firstVertexB: 11,
            faceA: 10, faceB: 11,
            prevA: 9, nextA: 28,
            prevB: 17, nextB: 20
        } // 29
    ];
    return [icosahedron, null];
}
function vectorAngle(first, second) {
    return Math.acos((first[0] * second[0] + first[1] * second[1] + first[2] * second[2]) / (Math.sqrt(first[0] * first[0] + first[1] * first[1] + first[2] * first[2]) * Math.sqrt(second[0] * second[0] + second[1] * second[1] + second[2] * second[2])));
}
function vectorLength(vector) {
    return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}
