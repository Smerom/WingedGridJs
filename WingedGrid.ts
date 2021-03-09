/* WARNING, CURRENTLY CAN ONLY SUBDIVIDE TRIANGULAR GRIDS. SUBDIVIDING DUALS IS UNDEFINED */

export class WingedEdge {
	grid: WingedGrid;
	index: number;

	firstVertexA(): number {
		return this.grid.edges_firstVertexA[this.index];
	}
	firstVertexB(): number {
		return this.grid.edges_firstVertexB[this.index];
	}
	faceA(): number {
		return this.grid.edges_faceA[this.index];
	}
	faceB(): number {
		return this.grid.edges_faceB[this.index];
	}
	prevA(): number {
		return this.grid.edges_prevA[this.index];
	}
	nextA(): number {
		return this.grid.edges_nextA[this.index];
	}
	prevB(): number {
		return this.grid.edges_prevB[this.index];
	}
	nextB(): number {
		return this.grid.edges_nextB[this.index];
	}

	constructor(sourceGrid: WingedGrid, edgeIndex: number) {
		this.grid = sourceGrid;
		this.index = edgeIndex;
	}
	/********* Base **********/
	nextEdgeForFace(faceIndex: number): [number, Error] {
		if (this.faceA() == faceIndex) {
			return [this.nextA(), null];
		} else if (this.faceB() == faceIndex) {
			return [this.nextB(), null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	prevEdgeForFace(faceIndex: number): [number, Error] {
		if (this.faceA() == faceIndex) {
			return [this.prevA(), null];
		} else if (this.faceB() == faceIndex) {
			return [this.prevB(), null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	firstVertexForFace(faceIndex: number): [number, Error] {
		if (this.faceA() == faceIndex) {
			return [this.firstVertexA(), null];
		} else if (this.faceB() == faceIndex) {
			return [this.firstVertexB(), null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	secondVertexForFace(faceIndex: number): [number, Error] {
		if (this.faceA() == faceIndex) {
			return [this.firstVertexB(), null];
		} else if (this.faceB() == faceIndex) {
			return [this.firstVertexA(), null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	nextEdgeForVertex(vertexIndex: number): [number, Error] {
		if (this.firstVertexA() == vertexIndex) {
			return [this.prevA(), null];
		} else if (this.firstVertexB() == vertexIndex) {
			return [this.prevB(), null];
		}
		return [-1, Error("Edge not associated with vertex.")];
	}
	prevEdgeForVertex(vertexIndex: number): [number, Error] {
		if (this.firstVertexA() == vertexIndex) {
			return this.nextEdgeForFace(this.faceB());
		} else if (this.firstVertexB() == vertexIndex) {
			return this.nextEdgeForFace(this.faceA());
		}
		return [-1, Error("Edge not associated with vertex.")];
	}
	adjacentForFace(faceIndex: number): [number, Error] {
		if (this.faceA() == faceIndex) {
			return [this.faceB(), null];
		} else if (this.faceB() == faceIndex) {
			return [this.faceA(), null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	adjacentForVertex(vertexIndex: number): [number, Error] {
		if (this.firstVertexA() == vertexIndex) {
			return [this.firstVertexB(), null]
		} else if (this.firstVertexB() == vertexIndex) {
			return [this.firstVertexA(), null]
		}
		return [-1, Error("Edge not associated with vertex.")];
	}
}

export class WingedFace {
	grid: WingedGrid;
	index: number;

	edges(): Array<WingedEdge> {
		let end: number;
		if (this.index == this.grid.faces_edgeOffsets.length - 1) {
			end = this.grid.faces_edges.length;
		} else {
			end = this.grid.faces_edgeOffsets[this.index + 1];
		}
		let edges = new Array<WingedEdge>(end - this.grid.faces_edgeOffsets[this.index]);
		let i: number = 0;
		for (var faceEdge = this.grid.faces_edgeOffsets[this.index]; faceEdge < end; ++faceEdge) {
			edges[i] = new WingedEdge(this.grid, this.grid.faces_edges[faceEdge]);
			i++;
		}
		return edges;
	}
	constructor(sourceGrid: WingedGrid, faceIndex: number) {
		this.grid = sourceGrid;
		this.index = faceIndex;
	}
}

export class WingedVertex {
	grid: WingedGrid;
	index: number;

	coords(): [number, number, number] {
		return [this.grid.vertices_coords[this.index * 3], this.grid.vertices_coords[this.index * 3 + 1], this.grid.vertices_coords[this.index * 3 + 2]];
	};
	edges(): Array<WingedEdge> {
		let end: number;
		if (this.index == this.grid.vertices_edgeOffsets.length - 1) {
			end = this.grid.vertices_edges.length;
		} else {
			end = this.grid.vertices_edgeOffsets[this.index + 1];
		}
		let edges = new Array<WingedEdge>(end - this.grid.vertices_edgeOffsets[this.index]);
		let i: number = 0;
		for (var faceEdge = this.grid.vertices_edgeOffsets[this.index]; faceEdge < end; ++faceEdge) {
			edges[i] = new WingedEdge(this.grid, this.grid.vertices_edges[faceEdge]);
			i++;
		}
		return edges;
	}
	vertexNeighbors(): Array<WingedVertex> {
		let neighborsIndecies = this.grid.neighborsForVertex(this.index);
		let neighbors = new Array<WingedVertex>(neighborsIndecies.length);
		for (var i = neighborsIndecies.length - 1; i >= 0; i--) {
			neighbors[i] = new WingedVertex(this.grid, neighborsIndecies[i]);
		}
		return neighbors;
	};

	constructor(sourceGrid: WingedGrid, faceIndex: number) {
		this.grid = sourceGrid;
		this.index = faceIndex;
	}
}

export class WingedGrid {
	faces(index: number): WingedFace {
		return new WingedFace(this, index);
	}
	faces_edgeOffsets: Int32Array;
	faces_edges: Int32Array;
	edges(index: number): WingedEdge {
		return new WingedEdge(this, index);
	}
	edges_firstVertexA: Int32Array;
	edges_firstVertexB: Int32Array;
	edges_faceA: Int32Array;
	edges_faceB: Int32Array;
	edges_prevA: Int32Array;
	edges_nextA: Int32Array;
	edges_prevB: Int32Array;
	edges_nextB: Int32Array;
	vertices(index: number): WingedVertex {
		return new WingedVertex(this, index);
	}
	vertices_coords: Float32Array;
	vertices_edgeOffsets: Int32Array;
	vertices_edges: Int32Array;
	constructor(faceCount: number, edgeCount: number, vertCount: number, totalFaceEdgeCount: number, totalVertexEdgeCount: number) {
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
	/********* Base **********/
	// TODO, OPTIMIZE FOR NEW STRUCTURE
	neighborsForFace(faceIndex: number): Array<number> {
		let neighbors = Array<number>();
		for (let i = this.faces(faceIndex).edges().length - 1; i >= 0; i--) {
			let [neighbor, err] = this.faces(faceIndex).edges()[i].adjacentForFace(faceIndex);
			neighbors[i] = neighbor;
		}
		return neighbors;
	}
	neighborsForVertex(vertexIndex: number): Array<number> {
		let neighbors = Array<number>();
		for (var i = this.vertices(vertexIndex).edges().length - 1; i >= 0; i--) {
			let [neighbor, error] = this.vertices(vertexIndex).edges()[i].adjacentForVertex(vertexIndex);
			neighbors[i] = neighbor;
		}
		return neighbors
	}

	/*********** SPHERE ************/
	normalizeVerticesToDistanceFromOrigin(wantedLength: number) {
		for (var i = this.vertices_edgeOffsets.length - 1; i >= 0; i--) {
			let vertex = this.vertices(i);
			let currentLength = vectorLength(vertex.coords());
			this.vertices_coords[i * 3 + 0] = this.vertices_coords[i * 3 + 0] * wantedLength / currentLength;
			this.vertices_coords[i * 3 + 1] = this.vertices_coords[i * 3 + 1] * wantedLength / currentLength;
			this.vertices_coords[i * 3 + 2] = this.vertices_coords[i * 3 + 2] * wantedLength / currentLength;
		}
	}

	/******** DUAL *********/
	createDual(): [WingedGrid, Error] {
		let dualGrid = new WingedGrid(this.vertices_edgeOffsets.length, this.edges_firstVertexA.length, this.faces_edgeOffsets.length, 0, 0); // swap vert and face counts
		// create faces, just swap with old vert edges
		dualGrid.faces_edges = this.vertices_edges;
		dualGrid.faces_edgeOffsets = this.vertices_edgeOffsets;


		// create vertices
		// swap with old face edges
		dualGrid.vertices_edges = this.faces_edges;
		dualGrid.vertices_edgeOffsets = this.faces_edgeOffsets;
		for (var i = this.faces_edgeOffsets.length - 1; i >= 0; i--) {
			let faceEdges = this.faces(i).edges();
			// set coords from center of old face
			let faceCenter: [number, number, number] = [0, 0, 0];
			let count: number = 0;
			for (var j = faceEdges.length - 1; j >= 0; j--) {
				let edge = faceEdges[j];
				let [vertexIndex, err] = edge.firstVertexForFace(i);
				if (err != null) {
					return [dualGrid, err]
				}
				faceCenter[0] += this.vertices_coords[vertexIndex * 3 + 0];
				faceCenter[1] += this.vertices_coords[vertexIndex * 3 + 1];
				faceCenter[2] += this.vertices_coords[vertexIndex * 3 + 2];
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
		for (let faceIndex = dualGrid.faces.length - 1; faceIndex >= 0; faceIndex--) {
			let faceEdges = dualGrid.faces(faceIndex).edges();
			for (var faceEdgeIndex = faceEdges.length - 1; faceEdgeIndex >= 0; faceEdgeIndex--) {
				let edgeIndex = faceEdges[faceEdgeIndex].index;
				if (dualGrid.edges_faceA[edgeIndex] == faceIndex) {
					if (faceEdgeIndex == 0) {
						dualGrid.edges_prevA[edgeIndex] = faceEdges[faceEdges.length - 1].index;
						dualGrid.edges_nextA[edgeIndex] = faceEdges[1].index;
					} else if (faceEdgeIndex == faceEdges.length - 1) {
						dualGrid.edges_prevA[edgeIndex] = faceEdges[faceEdges.length - 1].index;
						dualGrid.edges_nextA[edgeIndex] = faceEdges[0].index;
					} else {
						dualGrid.edges_prevA[edgeIndex] = faceEdges[faceEdgeIndex - 1].index;
						dualGrid.edges_nextA[edgeIndex] = faceEdges[faceEdgeIndex + 1].index;
					}
				}
				if (dualGrid.edges_faceB[edgeIndex] == faceIndex) {
					if (faceEdgeIndex == 0) {
						dualGrid.edges_prevB[edgeIndex] = faceEdges[faceEdges.length - 1].index;
						dualGrid.edges_nextB[edgeIndex] = faceEdges[1].index;
					} else if (faceEdgeIndex == faceEdges.length - 1) {
						dualGrid.edges_prevB[edgeIndex] = faceEdges[faceEdges.length - 1].index;
						dualGrid.edges_nextB[edgeIndex] = faceEdges[0].index;
					} else {
						dualGrid.edges_prevB[edgeIndex] = faceEdges[faceEdgeIndex - 1].index;
						dualGrid.edges_nextB[edgeIndex] = faceEdges[faceEdgeIndex + 1].index;
					}
				}
			}
		}
		return [dualGrid, null]
	}

	/************* Subdivision helpers ****************/
	// rename vertexInexAtClockwiseIndexOnFaceToSubdivide
	vertexIndexAtClockwiseIndexOnOldFace(faceIndex: number, edgeInFaceIndex: number, clockwiseVertexIndex: number, edgeSubdivisions: number): number {
		let edge = this.faces(faceIndex).edges()[edgeInFaceIndex];
		let edgeIndex = edge.index;
		if (this.edges_faceA[edgeIndex] == faceIndex) {
			return this.vertices_edgeOffsets.length + edgeIndex * edgeSubdivisions + clockwiseVertexIndex;
		}
		if (this.edges_faceB[edgeIndex] == faceIndex) {
			return this.vertices_edgeOffsets.length + edgeIndex * edgeSubdivisions + edgeSubdivisions - 1 - clockwiseVertexIndex;
		}
		return -1;
	}
	edgeIndexAtClockwiseIndexOnOldFace(faceIndex: number, edgeInFaceIndex: number, clockwiseEdgeIndex: number, edgeSubdivisions: number): number {
		let oldEdge = this.faces(faceIndex).edges()[edgeInFaceIndex];
		let oldEdgeIndex = oldEdge.index;
		if (this.edges_faceA[oldEdgeIndex] == faceIndex) {
			return oldEdgeIndex * (edgeSubdivisions + 1) + clockwiseEdgeIndex;
		}
		if (this.edges_faceB[oldEdgeIndex] == faceIndex) {
			return oldEdgeIndex * (edgeSubdivisions + 1) + edgeSubdivisions - clockwiseEdgeIndex;
		}
		return -1;
	}

	/*************** SUBDIVISIONS ***************/
	subdivideTriangles(edgeSubdivisions: number): [WingedGrid, Error] {
		if (edgeSubdivisions < 1) {
			return [null, Error("Invalid number of subdivisions.")]
		}
		// subdividing each edge n times produces a number of faces equal
		//  to the base face count multiplied by (1/2(n+2)(n+1) + 1/2(n+1)(n))
		const faceCount: number = this.faces_edgeOffsets.length * ((edgeSubdivisions + 2) * (edgeSubdivisions + 1) / 2 + (edgeSubdivisions + 1) * edgeSubdivisions / 2);
		let dividedGrid = new WingedGrid(
			faceCount,
			// since each face 'owns' 1/2 of three edges, there are 1.5 times as
			//  many edges as faces
			3 * faceCount / 2,
			// Euler-somebody or other gives us the vertex count of
			//  1/2*(face count) + 2
			faceCount / 2 + 2,
			faceCount * 3, // three for each face
			(faceCount / 2 + 2 - this.vertices_edgeOffsets.length) * 6 + this.vertices_edges.length
		);

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
		let currentOffset: number = 0;
		for (var i = 0; i < dividedGrid.vertices_edgeOffsets.length; i++) {
			if (i < this.vertices_edgeOffsets.length) {
				currentOffset = this.vertices_edgeOffsets[i];
			} else if (i == this.vertices_edgeOffsets.length) {
				currentOffset = this.vertices_edges.length;
			} else {
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
		let origVertexCount = this.vertices_edgeOffsets.length;
		let origEdgeCount = this.edges_faceA.length;
		for (let i = origEdgeCount - 1; i >= 0; i--) {
			// first edge has origional vertex
			dividedGrid.edges_firstVertexA[i * (edgeSubdivisions + 1)] = this.edges_firstVertexA[i];
			dividedGrid.edges_firstVertexB[i * (edgeSubdivisions + 1)] = origVertexCount + i * edgeSubdivisions;

			for (let j = 1; j < edgeSubdivisions; j++) {
				// set the edge vertex indecies
				dividedGrid.edges_firstVertexA[i * (edgeSubdivisions + 1) + j] = origVertexCount + i * edgeSubdivisions + j - 1;
				dividedGrid.edges_firstVertexB[i * (edgeSubdivisions + 1) + j] = origVertexCount + i * edgeSubdivisions + j;
			}

			// connect last new edge
			dividedGrid.edges_firstVertexA[i * (edgeSubdivisions + 1) + edgeSubdivisions] = origVertexCount + i * edgeSubdivisions + edgeSubdivisions - 1;
			dividedGrid.edges_firstVertexB[i * (edgeSubdivisions + 1) + edgeSubdivisions] = this.edges_firstVertexB[i];
		}

		/********* Edges created interior to old faces go in the second section, ordered by face. ****/
		for (let faceIndex = this.faces_edgeOffsets.length - 1; faceIndex >= 0; faceIndex--) {
			let oldFace = this.faces(faceIndex);
			// vertex offset for vertices interior to the face
			let vertexOffset = origVertexCount + origEdgeCount * edgeSubdivisions + (edgeSubdivisions * (edgeSubdivisions - 1) / 2) * faceIndex;

			var edgeOffset: number = (edgeSubdivisions + 1) * origEdgeCount + 3 * (edgeSubdivisions) * (edgeSubdivisions + 1) / 2 * faceIndex;
			if (edgeSubdivisions == 1) {
				dividedGrid.edges_firstVertexA[edgeOffset + 0] = origVertexCount + oldFace.edges()[0].index;
				dividedGrid.edges_firstVertexB[edgeOffset + 0] = origVertexCount + oldFace.edges()[1].index;

				dividedGrid.edges_firstVertexA[edgeOffset + 1] = origVertexCount + oldFace.edges()[0].index;
				dividedGrid.edges_firstVertexB[edgeOffset + 1] = origVertexCount + oldFace.edges()[2].index;

				dividedGrid.edges_firstVertexA[edgeOffset + 2] = origVertexCount + oldFace.edges()[1].index;
				dividedGrid.edges_firstVertexB[edgeOffset + 2] = origVertexCount + oldFace.edges()[2].index;
			} else {
				// first row
				dividedGrid.edges_firstVertexA[edgeOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions - 1, edgeSubdivisions);
				dividedGrid.edges_firstVertexB[edgeOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 0, edgeSubdivisions);

				dividedGrid.edges_firstVertexA[edgeOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions - 1, edgeSubdivisions);
				dividedGrid.edges_firstVertexB[edgeOffset + 1] = vertexOffset + 0;

				dividedGrid.edges_firstVertexA[edgeOffset + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 0, edgeSubdivisions);
				dividedGrid.edges_firstVertexB[edgeOffset + 2] = vertexOffset + 0;

				// middle rows
				var rowOffset: number;
				var i: number = 0;
				for (i = 1; i < edgeSubdivisions - 1; i++) {
					rowOffset = i * (i + 1) * 3 / 2
					// first border
					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions - 1 - i, edgeSubdivisions);
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 0] = vertexOffset + (i * (i - 1) / 2);

					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions - 1 - i, edgeSubdivisions);
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 1] = vertexOffset + (i * (i - 1) / 2) + i;

					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 2] = vertexOffset + (i * (i - 1) / 2);
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 2] = vertexOffset + (i * (i - 1)) / 2 + i;
					// interior of face
					for (let j = 1; j < i; j++) {
						dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset + (i * (i - 1) / 2) + j - 1;
						dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset + (i * (i - 1) / 2) + j;

						dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 1] = vertexOffset + (i * (i - 1) / 2) + j - 1;
						dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 1] = vertexOffset + (i * (i - 1) / 2) + i + j;

						dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 2] = vertexOffset + (i * (i - 1) / 2) + j;
						dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 2] = vertexOffset + (i * (i - 1) / 2) + i + j;
					}

					// second border
					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 0] = vertexOffset + (i * (i - 1) / 2) + i - 1;
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);

					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 1] = vertexOffset + (i * (i - 1) / 2) + i - 1;
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 1] = vertexOffset + (i * (i - 1) / 2) + 2 * i;

					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 2] = vertexOffset + (i * (i - 1) / 2) + 2 * i;
				}

				// last row
				i = edgeSubdivisions - 1; // should already be set, but just incase
				rowOffset = i * (i + 1) * 3 / 2;
				// border 1
				dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, 0, edgeSubdivisions);
				dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 0] = vertexOffset + (i * (i - 1) / 2);

				dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, 0, edgeSubdivisions);
				dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions - 1, edgeSubdivisions);

				dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + 2] = vertexOffset + (i * (i - 1) / 2);
				dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions - 1, edgeSubdivisions);

				// middle
				for (let j = 1; j < i; j++) {
					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset + (i * (i - 1) / 2) + j - 1;
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 0] = vertexOffset + (i * (i - 1) / 2) + j;

					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 1] = vertexOffset + (i * (i - 1) / 2) + j - 1;
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions - j - 1, edgeSubdivisions);

					dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + j * 3 + 2] = vertexOffset + (i * (i - 1) / 2) + j;
					dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + j * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions - j - 1, edgeSubdivisions);
				}

				// border 2
				dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 0] = vertexOffset + (i * (i - 1) / 2) + i - 1;
				dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 0] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);

				dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 1] = vertexOffset + (i * (i - 1) / 2) + i - 1;
				dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 1] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, 0, edgeSubdivisions);

				dividedGrid.edges_firstVertexA[edgeOffset + rowOffset + i * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);
				dividedGrid.edges_firstVertexB[edgeOffset + rowOffset + i * 3 + 2] = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, 0, edgeSubdivisions);
			}
		}

		// create the faces, update the edges
		/******************* FACE SUBDIVISION ***********************/
		for (let faceIndex = this.faces_edgeOffsets.length - 1; faceIndex >= 0; faceIndex--) {
			// get the number of faces we divide this one into
			let subFaceCount = (edgeSubdivisions + 2) * (edgeSubdivisions + 1) / 2 + (edgeSubdivisions + 1) * edgeSubdivisions / 2;

			// number of internal edges created on this face
			//  equal to ( (edgeSubdivisions + 1) choose 2 ) * 3
			var subEdgeCount = 3 * edgeSubdivisions * (edgeSubdivisions + 1) / 2;

			// each new face set starts at
			var indexStart: number = subFaceCount * faceIndex;

			var edgeOffset: number = (edgeSubdivisions + 1) * origEdgeCount + subEdgeCount * faceIndex;

			// first corner
			// set edges in face
			dividedGrid.faces_edgeOffsets[indexStart] = (indexStart) * 3; // three edges per face
			dividedGrid.faces_edges[indexStart * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions, edgeSubdivisions);
			dividedGrid.faces_edges[indexStart * 3 + 1] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 0, edgeSubdivisions);
			dividedGrid.faces_edges[indexStart * 3 + 2] = edgeOffset;

			// loop through the middle section of faces
			// edges grow by 1/2*i*(i-1)*3
			for (let i = 1; i < edgeSubdivisions; i++) {
				var rowIndexStart: number = i * (i + 1) / 2 + i * (i - 1) / 2

				// edge
				dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart] = (indexStart + rowIndexStart) * 3;
				dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions - i, edgeSubdivisions);
				dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 1] = edgeOffset + i * (i - 1) * 3 / 2 + 1;
				dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 2] = edgeOffset + i * (i + 1) * 3 / 2;

				// middle
				// up to one less than next index start
				for (let j = 1; j < (i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - rowIndexStart - 1; j++) {
					if (j % 2 == 1) {
						dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
						dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = edgeOffset + i * (i - 1) * 3 / 2 + (j - 1) * 3 / 2;
						dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + i * (i - 1) * 3 / 2 + (j - 1) * 3 / 2 + 2;
						dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + i * (i - 1) * 3 / 2 + (j - 1) * 3 / 2 + 1;
					} else {
						dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3
						dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = edgeOffset + i * (i + 1) * 3 / 2 + j * 3 / 2;
						dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + i * (i - 1) * 3 / 2 + (j - 2) * 3 / 2 + 2;
						dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + i * (i - 1) * 3 / 2 + j * 3 / 2 + 1;
					}
				}

				// edge
				dividedGrid.faces_edgeOffsets[indexStart + (i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - 1] = (indexStart + (i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - 1) * 3;
				dividedGrid.faces_edges[(indexStart + (i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - 1) * 3 + 0] = edgeOffset + i * (i + 1) * 3 / 2 + ((i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - rowIndexStart - 1) * 3 / 2;
				dividedGrid.faces_edges[(indexStart + (i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - 1) * 3 + 1] = edgeOffset + i * (i - 1) * 3 / 2 + ((i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - rowIndexStart - 3) * 3 / 2 + 2;
				dividedGrid.faces_edges[(indexStart + (i + 1) * (i + 2) / 2 + (i + 1) * i / 2 - 1) * 3 + 2] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);
			}

			var rowIndexStart = edgeSubdivisions * (edgeSubdivisions + 1) / 2 + edgeSubdivisions * (edgeSubdivisions - 1) / 2
			// bottom corner 1
			dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart] = (indexStart + rowIndexStart) * 3;
			dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 0, 0, edgeSubdivisions);
			dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + 1;
			dividedGrid.faces_edges[(indexStart + rowIndexStart) * 3 + 2] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions, edgeSubdivisions);
			// bottom edge
			for (let j = 1; j < subFaceCount - rowIndexStart - 1; j++) {
				if (j % 2 == 1) {
					dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
					dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2;
					dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2 + 2;
					dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 1) * 3 / 2 + 1;
				} else {
					dividedGrid.faces_edgeOffsets[indexStart + rowIndexStart + j] = (indexStart + rowIndexStart + j) * 3;
					dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions - j / 2, edgeSubdivisions);
					dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (j - 2) * 3 / 2 + 2;
					dividedGrid.faces_edges[(indexStart + rowIndexStart + j) * 3 + 2] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + j * 3 / 2 + 1;
				}
			}
			// bottom corner 2
			dividedGrid.faces_edgeOffsets[indexStart + subFaceCount - 1] = (indexStart + subFaceCount - 1) * 3;
			dividedGrid.faces_edges[(indexStart + subFaceCount - 1) * 3 + 0] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 2, 0, edgeSubdivisions);
			dividedGrid.faces_edges[(indexStart + subFaceCount - 1) * 3 + 1] = edgeOffset + edgeSubdivisions * (edgeSubdivisions - 1) * 3 / 2 + (subFaceCount - rowIndexStart - 3) * 3 / 2 + 2;
			dividedGrid.faces_edges[(indexStart + subFaceCount - 1) * 3 + 2] = this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 1, edgeSubdivisions, edgeSubdivisions);
		}

		// set edge faces from the previously build edge arrays
		for (var faceIndex = dividedGrid.faces_edgeOffsets.length - 1; faceIndex >= 0; faceIndex--) {
			let edges = dividedGrid.faces(faceIndex).edges();
			let edgesLength = edges.length;
			let thisEdge: number, nextEdge: number, prevEdge: number;

			prevEdge = edges[edgesLength - 1].index;
			thisEdge = edges[0].index;
			nextEdge = edges[1].index;
			// test vertices
			if (dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexA[prevEdge] || dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexB[prevEdge]) {
				// check the next edge also matches and face A is not set
				if (dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexA[nextEdge] || dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexB[nextEdge]) {
					if (dividedGrid.edges_faceA[thisEdge] == -1) {
						dividedGrid.edges_faceA[thisEdge] = faceIndex;
						dividedGrid.edges_prevA[thisEdge] = prevEdge;
						dividedGrid.edges_nextA[thisEdge] = nextEdge;
					} else {
						console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge + " With edge set: " + edges);
					}
				} else {
					console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else if (dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexA[prevEdge] || dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexB[prevEdge]) {
				// check the next edge also matches and face B is not set
				if (dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexA[nextEdge] || dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexB[nextEdge]) {
					if (dividedGrid.edges_faceB[thisEdge] == -1) {
						dividedGrid.edges_faceB[thisEdge] = faceIndex;
						dividedGrid.edges_prevB[thisEdge] = prevEdge;
						dividedGrid.edges_nextB[thisEdge] = nextEdge;
					} else {
						console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge + " With edge set: " + edges);
					}
				} else {
					console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else {
				console.log("For face " + faceIndex + ". Edges Don't share a vertex!");
			}

			// loop through middle edges
			for (let i = 1; i < edgesLength - 1; i++) {
				prevEdge = edges[i - 1].index;
				thisEdge = edges[i].index;
				nextEdge = edges[i + 1].index;
				// test vertcies
				if (dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexA[prevEdge] || dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexB[prevEdge]) {
					// check the next edge also matches and face A is not set
					if (dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexA[nextEdge] || dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexB[nextEdge]) {
						if (dividedGrid.edges_faceA[thisEdge] == -1) {
							dividedGrid.edges_faceA[thisEdge] = faceIndex;
							dividedGrid.edges_prevA[thisEdge] = prevEdge;
							dividedGrid.edges_nextA[thisEdge] = nextEdge;
						} else {
							console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge + " With edge set: " + edges);
						}
					} else {
						console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
					}
				} else if (dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexA[prevEdge] || dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexB[prevEdge]) {
					// check the next edge also matches and face B is not set
					if (dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexA[nextEdge] || dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexB[nextEdge]) {
						if (dividedGrid.edges_faceB[thisEdge] == -1) {
							dividedGrid.edges_faceB[thisEdge] = faceIndex;
							dividedGrid.edges_prevB[thisEdge] = prevEdge;
							dividedGrid.edges_nextB[thisEdge] = nextEdge;
						} else {
							console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge + " With edge set: " + edges);
						}
					} else {
						console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
					}
				} else {
					console.log("For face " + faceIndex + ". Edges Don't share a vertex!");
				}
			}

			// last edge
			prevEdge = edges[edgesLength - 2].index;
			thisEdge = edges[edgesLength - 1].index;
			nextEdge = edges[0].index;
			// test vertices
			if (dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexA[prevEdge] || dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexB[prevEdge]) {
				// check the next edge also matches and face A is not set
				if (dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexA[nextEdge] || dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexB[nextEdge]) {
					if (dividedGrid.edges_faceA[thisEdge] == -1) {
						dividedGrid.edges_faceA[thisEdge] = faceIndex;
						dividedGrid.edges_prevA[thisEdge] = prevEdge;
						dividedGrid.edges_nextA[thisEdge] = nextEdge;
					} else {
						console.log("For face " + faceIndex + ". Face A has already been set for edge: " + thisEdge + " With edge set: " + edges);
					}
				} else {
					console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else if (dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexA[prevEdge] || dividedGrid.edges_firstVertexB[thisEdge] == dividedGrid.edges_firstVertexB[prevEdge]) {
				// check the next edge also matches and face B is not set
				if (dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexA[nextEdge] || dividedGrid.edges_firstVertexA[thisEdge] == dividedGrid.edges_firstVertexB[nextEdge]) {
					if (dividedGrid.edges_faceB[thisEdge] == -1) {
						dividedGrid.edges_faceB[thisEdge] = faceIndex;
						dividedGrid.edges_prevB[thisEdge] = prevEdge;
						dividedGrid.edges_nextB[thisEdge] = nextEdge;
					} else {
						console.log("For face " + faceIndex + ". Face B has already been set for edge: " + thisEdge + " With edge set: " + edges);
					}
				} else {
					console.log("For face " + faceIndex + ". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else {
				console.log("For face " + faceIndex + ". Edges Don't share a vertex!");
			}//*/
		}

		// create the vertices
		/******************* VERTEX SUBDIVISION ***********************/
		// set coords for the origional verts
		// maybe copy slice directly?
		for (let i = this.vertices_edgeOffsets.length - 1; i >= 0; i--) {
			dividedGrid.vertices_coords[i * 3 + 0] = this.vertices_coords[i * 3 + 0];
			dividedGrid.vertices_coords[i * 3 + 1] = this.vertices_coords[i * 3 + 1];
			dividedGrid.vertices_coords[i * 3 + 2] = this.vertices_coords[i * 3 + 2];
		}

		// subdivide along each edge
		for (let i = this.edges_faceA.length - 1; i >= 0; i--) {
			let firstVertex = dividedGrid.vertices(this.edges_firstVertexA[i]);
			let firstVertexIndex = firstVertex.index;
			let secondVertex = dividedGrid.vertices(this.edges_firstVertexB[i]);
			let secondVertexIndex = secondVertex.index;
			// angle to subdivide
			let angleToSubdivide: number;
			angleToSubdivide = vectorAngle(firstVertex.coords(), secondVertex.coords())

			// angle between origin, first vertex, and second vertex
			let vectorA: [number, number, number] = [-1, -1, -1], vectorB: [number, number, number] = [-1, -1, -1];
			vectorA[0] = -1 * dividedGrid.vertices_coords[firstVertexIndex * 3 + 0];
			vectorA[1] = -1 * dividedGrid.vertices_coords[firstVertexIndex * 3 + 1];
			vectorA[2] = -1 * dividedGrid.vertices_coords[firstVertexIndex * 3 + 2];

			vectorB[0] = dividedGrid.vertices_coords[secondVertexIndex * 3 + 0] + vectorA[0];
			vectorB[1] = dividedGrid.vertices_coords[secondVertexIndex * 3 + 1] + vectorA[1];
			vectorB[2] = dividedGrid.vertices_coords[secondVertexIndex * 3 + 2] + vectorA[2];

			let cornerAngle: number;
			cornerAngle = vectorAngle(vectorA, vectorB);

			// unit vector from first to second vertex
			let stepDirection: [number, number, number] = [-1, -1, -1];
			stepDirection[0] = vectorB[0] / vectorLength(vectorB);
			stepDirection[1] = vectorB[1] / vectorLength(vectorB);
			stepDirection[2] = vectorB[2] / vectorLength(vectorB);

			// origional radius of the
			let sphereRadius: number = vectorLength(firstVertex.coords());

			let divisionLength: number;
			for (let j = 0; j < edgeSubdivisions; j++) {
				// find the new vertex position and create the vertex
				// but don't correct it's length yet
				divisionLength = Math.sin(angleToSubdivide * ((j + 1) / (edgeSubdivisions + 1))) * sphereRadius / Math.sin(Math.PI - cornerAngle - angleToSubdivide * ((j + 1) / (edgeSubdivisions + 1)));

				dividedGrid.vertices_coords[(origVertexCount + i * edgeSubdivisions + j) * 3 + 0] = dividedGrid.vertices_coords[firstVertexIndex * 3 + 0] + stepDirection[0] * divisionLength;
				dividedGrid.vertices_coords[(origVertexCount + i * edgeSubdivisions + j) * 3 + 1] = dividedGrid.vertices_coords[firstVertexIndex * 3 + 1] + stepDirection[1] * divisionLength;
				dividedGrid.vertices_coords[(origVertexCount + i * edgeSubdivisions + j) * 3 + 2] = dividedGrid.vertices_coords[firstVertexIndex * 3 + 2] + stepDirection[2] * divisionLength;
			}
		}

		// subdivide face interior
		for (let faceIndex = 0; faceIndex < this.faces_edgeOffsets.length; faceIndex++) {
			var vertexOffset: number = origVertexCount + origEdgeCount * edgeSubdivisions + edgeSubdivisions * (edgeSubdivisions - 1) / 2 * faceIndex;
			// only if we have more than one division
			if (edgeSubdivisions > 1) {
				for (let i = 0; i < edgeSubdivisions - 1; i++) {
					let firstVertex = dividedGrid.vertices(this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions - 2 - i, edgeSubdivisions));
					let firstVertexIndex = firstVertex.index;
					let secondVertex = dividedGrid.vertices(this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 1 + i, edgeSubdivisions));
					let secondVertexIndex = secondVertex.index;
					// angle to subdivide
					var angleToSubdivide: number;
					angleToSubdivide = vectorAngle(firstVertex.coords(), secondVertex.coords());

					// angle between origin, first vertex, and second vertex
					var vectorA: [number, number, number] = [-1, -1, -1], vectorB: [number, number, number] = [-1, -1, -1];
					vectorA[0] = -1 * dividedGrid.vertices_coords[firstVertexIndex * 3 + 0];
					vectorA[1] = -1 * dividedGrid.vertices_coords[firstVertexIndex * 3 + 1];
					vectorA[2] = -1 * dividedGrid.vertices_coords[firstVertexIndex * 3 + 2];

					vectorB[0] = dividedGrid.vertices_coords[secondVertexIndex * 3 + 0] + vectorA[0];
					vectorB[1] = dividedGrid.vertices_coords[secondVertexIndex * 3 + 1] + vectorA[1];
					vectorB[2] = dividedGrid.vertices_coords[secondVertexIndex * 3 + 2] + vectorA[2];

					var cornerAngle: number;
					cornerAngle = vectorAngle(vectorA, vectorB);

					// unit vector from first to second vertex
					var stepDirection: [number, number, number] = [-1, -1, -1];
					stepDirection[0] = vectorB[0] / vectorLength(vectorB);
					stepDirection[1] = vectorB[1] / vectorLength(vectorB);
					stepDirection[2] = vectorB[2] / vectorLength(vectorB);

					var firstVectorLength: number = vectorLength(firstVertex.coords());
					var divisionLength: number = 0
					for (let j = 0; j < i + 1; j++) {
						divisionLength = Math.sin(angleToSubdivide * ((j + 1) / (i + 2))) * firstVectorLength / Math.sin(Math.PI - cornerAngle - angleToSubdivide * ((j + 1) / (i + 2)));
						if (vertexOffset + (i * (i + 1) / 2) + j >= dividedGrid.vertices_edgeOffsets.length) {
							console.log("breakpoint");
						}
						dividedGrid.vertices_coords[(vertexOffset + (i * (i + 1) / 2) + j) * 3 + 0] = dividedGrid.vertices_coords[firstVertexIndex * 3 + 0] + stepDirection[0] * divisionLength;
						dividedGrid.vertices_coords[(vertexOffset + (i * (i + 1) / 2) + j) * 3 + 1] = dividedGrid.vertices_coords[firstVertexIndex * 3 + 1] + stepDirection[1] * divisionLength;
						dividedGrid.vertices_coords[(vertexOffset + (i * (i + 1) / 2) + j) * 3 + 2] = dividedGrid.vertices_coords[firstVertexIndex * 3 + 2] + stepDirection[2] * divisionLength;
					}
				}
			}
		}

		// set vertex edge array
		// loop through edges so we only have to touch each one once
		for (let index = dividedGrid.edges_faceA.length - 1; index >= 0; index--) {
			var nextEdgeIndex: number = -1;
			var theVertexIndex: number = dividedGrid.edges_firstVertexA[index];
			var error: Error;
			if (dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[0]] == -1) {
				if (dividedGrid.edges_firstVertexA[index] == theVertexIndex) {
					nextEdgeIndex = dividedGrid.edges_prevA[index];
				} else if (dividedGrid.edges_firstVertexB[index] == theVertexIndex) {
					nextEdgeIndex = dividedGrid.edges_prevB[index];
				}
				dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertexIndex] + 0] = index;
				var i: number = 1;
				while (index != nextEdgeIndex) {
					dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertexIndex] + i] = nextEdgeIndex;
					if (dividedGrid.edges_firstVertexA[nextEdgeIndex] == theVertexIndex) {
						nextEdgeIndex = dividedGrid.edges_prevA[nextEdgeIndex];
					} else if (dividedGrid.edges_firstVertexB[nextEdgeIndex] == theVertexIndex) {
						nextEdgeIndex = dividedGrid.edges_prevB[nextEdgeIndex];
					}
					i = i + 1;
				}
			}
			theVertexIndex = dividedGrid.edges_firstVertexB[index];
			if (dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[0]] == -1) {
				if (dividedGrid.edges_firstVertexA[index] == theVertexIndex) {
					nextEdgeIndex = dividedGrid.edges_prevA[index];
				} else if (dividedGrid.edges_firstVertexB[index] == theVertexIndex) {
					nextEdgeIndex = dividedGrid.edges_prevB[index];
				}
				dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertexIndex] + 0] = index
				var i: number = 1;
				while (index != nextEdgeIndex) {
					dividedGrid.vertices_edges[dividedGrid.vertices_edgeOffsets[theVertexIndex] + i] = nextEdgeIndex;
					if (dividedGrid.edges_firstVertexA[nextEdgeIndex] == theVertexIndex) {
						nextEdgeIndex = dividedGrid.edges_prevA[nextEdgeIndex];
					} else if (dividedGrid.edges_firstVertexB[nextEdgeIndex] == theVertexIndex) {
						nextEdgeIndex = dividedGrid.edges_prevB[nextEdgeIndex];
					}
					i = i + 1;
				}
			}
		}

		return [dividedGrid, null];
	}
}

export function baseIcosahedron(): WingedGrid {
	let icosahedron = new WingedGrid(0, 0, 0, 0, 0); // arrays will be set manually

	icosahedron.faces_edges = Int32Array.from([
		// cap 1
		9, 4, 0 // 0
		, 1, 5, 0 // 1
		, 2, 6, 1 // 2
		, 3, 7, 2 // 3
		, 4, 8, 3 // 4
		// cap 2
		, 14, 19, 10 // 5
		, 15, 11, 10 // 6
		, 16, 12, 11 // 7
		, 17, 13, 12 // 8
		, 18, 14, 13 // 9
		// Ring of 10 between caps, clockwise around cap 1
		//   starting with the first face adjacent to face 0
		, 29, 28, 9 // 10
		, 29, 20, 17 // 11
		, 21, 20, 5 // 12
		, 21, 22, 18 // 13
		, 23, 22, 6 // 14
		, 23, 24, 19 // 15
		, 25, 24, 7 // 16
		, 25, 26, 15 // 17
		, 27, 26, 8 // 18
		, 27, 28, 16 // 19
	]);

	// set offsets, three apiece
	icosahedron.faces_edgeOffsets = new Int32Array(20);
	for (var i = 0; i < 20; ++i) {
		icosahedron.faces_edgeOffsets[i] = i * 3;
	}

	const goldenRatio = 1.61803398875
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

	return icosahedron
}

function vectorAngle(first: [number, number, number], second: [number, number, number]): number {
	return Math.acos((first[0] * second[0] + first[1] * second[1] + first[2] * second[2]) / (Math.sqrt(first[0] * first[0] + first[1] * first[1] + first[2] * first[2]) * Math.sqrt(second[0] * second[0] + second[1] * second[1] + second[2] * second[2])));
}

function vectorLength(vector: [number, number, number]): number {
	return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}



type Vector = [number, number, number]

export class CalculatedGrid {
	baseGrid: WingedGrid
	baseVerts: number
	baseEdges: number
	baseFaces: number

	constructor(base: WingedGrid) {
		this.baseGrid = base
		this.baseVerts = this.baseGrid.vertices_edgeOffsets.length
		this.baseEdges = this.baseGrid.edges_faceA.length
		this.baseFaces = this.baseGrid.faces_edgeOffsets.length
	}
	VertexAndNeighbors(idx: number, subDivs: number): [Vector, number[]] {

		const neighbors = this.VertexNeighbors(idx, subDivs)
		const vertex = this.Vertex(idx, subDivs)

		return [vertex, neighbors]
	}

	VertexNeighbors(idx: number, subDivs: number): number[] {
		const baseVerts = this.baseVerts
		const baseEdges = this.baseEdges
		const baseFaces = this.baseFaces
		// three cases
		// origional vert
		// vert on origional edge
		// vert inside original face
		if (idx < baseVerts) {
			// origional vert
			return this.vertexNeighborsOrigional(idx, subDivs)
		}
		idx -= baseVerts
		if (idx < baseEdges * subDivs) {
			return this.vertexNeighborsEdge(idx, subDivs)
		}
		idx -= baseEdges * subDivs
		if (idx < (subDivs - 1) * subDivs / 2 * baseFaces) {
			return this.vertexNeighborsFace(idx, subDivs)
		}
		throw ("past end")
	}

	vertexNeighborsOrigional(idx: number, subDivs: number): number[] {
		let neighbors = Array<number>(this.baseGrid.vertices(idx).edges.length)
		let origVertexCount = this.baseVerts

		this.baseGrid.vertices(idx).edges().forEach((edge, i) => {
			const edx = edge.index

			if (edge.firstVertexA() === idx) {
				neighbors[i] = origVertexCount + edx * subDivs
			} else {
				neighbors[i] = origVertexCount + edx * subDivs + subDivs - 1
			}
		});

		return neighbors
	}

	vertexNeighborsEdge(idx: number, subDivs: number): number[] {
		const baseVerts = this.baseVerts
		// six for all but 12 origional pentagons
		let neighbors = Array<number>(6)

		const [edgeIdx, div] = divmod(idx, subDivs)
		const edge = this.baseGrid.edges(edgeIdx)

		// two cases, right next to corner, and not
		if (div == 0 || div == subDivs - 1) {
			// on edge ones
			if (div == 0) {
				neighbors[0] = edge.firstVertexA()
				neighbors[3] = idx + baseVerts + 1
			} else {
				neighbors[0] = edge.firstVertexB()
				neighbors[3] = idx + baseVerts - 1
			}

			// find this edge on vertex
			const edges = this.baseGrid.vertices(neighbors[0]).edges()
			let edgeInVertEdx = -1
			for (let i = 0; i < edges.length; i++) {
				const edx = edges[i].index;
				if (edx == edgeIdx) {
					edgeInVertEdx = i
					break
				}
			}

			// 5 is next, 1 is previous
			const nextIdx = edges[(edges.length + edgeInVertEdx + 1) % edges.length].index
			const nextEdge = this.baseGrid.edges(nextIdx)
			if (nextEdge.firstVertexA() == neighbors[0]) {
				neighbors[5] = baseVerts + nextIdx * subDivs
			} else {
				neighbors[5] = baseVerts + nextIdx * subDivs + subDivs - 1
			}
			const prevIdx = edges[(edges.length + edgeInVertEdx - 1) % edges.length].index
			const prevEdge = this.baseGrid.edges(prevIdx)
			if (prevEdge.firstVertexA() == neighbors[0]) {
				neighbors[1] = baseVerts + prevIdx * subDivs
			} else {
				neighbors[1] = baseVerts + prevIdx * subDivs + subDivs - 1
			}

			// 2 and 4
			// offset from face options are
			// 0, (subdivs-2)*(subdivs-1), subdivs*(subdivs-1) - 1

			// check face a
			let faceFourIdx: number, faceTwoIdx: number
			if (div == 0) {
				faceFourIdx = edge.faceA()
				faceTwoIdx = edge.faceB()
			} else {
				faceFourIdx = edge.faceB()
				faceTwoIdx = edge.faceA()
			}
			const faceFour = this.baseGrid.faces(faceFourIdx)
			let offset: number
			if (faceFour.edges()[0].index == edgeIdx) {
				// we are next to first vertex, first element of last row
				offset = (subDivs - 2) * (subDivs - 1) / 2
			} else if (faceFour.edges()[1].index == edgeIdx) {
				// we are next to second vertex, first element
				offset = 0
			} else {
				// we are next to third vertext, last element
				offset = subDivs * (subDivs - 1) / 2 - 1
			}

			neighbors[4] = baseVerts + this.baseEdges * subDivs + (subDivs - 1) * subDivs / 2 * faceFourIdx + offset

			// check face a
			const faceTwo = this.baseGrid.faces(faceTwoIdx)
			if (faceTwo.edges()[0].index == edgeIdx) {
				// we are next to second vertex, first element
				offset = 0
			} else if (faceTwo.edges()[1].index == edgeIdx) {
				// we are next to third vertext, last element
				offset = subDivs * (subDivs - 1) / 2 - 1
			} else {
				// we are next to first vertex, first element of last row
				offset = (subDivs - 2) * (subDivs - 1) / 2
			}
			neighbors[2] = baseVerts + this.baseEdges * subDivs + (subDivs - 1) * subDivs / 2 * faceTwoIdx + offset

		} else {
			// not next to corner
			neighbors[0] = idx + baseVerts - 1
			neighbors[3] = idx + baseVerts + 1

			const faceA = this.baseGrid.faces(edge.faceA())
			let faceOffset = baseVerts + this.baseEdges * subDivs + (subDivs - 1) * subDivs / 2 * edge.faceA()

			if (faceA.edges()[0].index == edgeIdx) {
				neighbors[4] = faceOffset + firstInRow(subDivs - div - 1)
				neighbors[5] = faceOffset + firstInRow(subDivs - div)
			} else if (faceA.edges()[1].index == edgeIdx) {
				neighbors[4] = faceOffset + lastInRow(div + 1)
				neighbors[5] = faceOffset + lastInRow(div)
			} else {
				// parallel to rows
				neighbors[4] = faceOffset + lastInRow(subDivs - 1) - div
				neighbors[5] = faceOffset + lastInRow(subDivs - 1) - div + 1
			}

			const faceB = this.baseGrid.faces(edge.faceB())
			faceOffset = baseVerts + this.baseEdges * subDivs + (subDivs - 1) * subDivs / 2 * edge.faceB()
			if (faceB.edges()[0].index == edgeIdx) {
				neighbors[1] = faceOffset + firstInRow(div)
				neighbors[2] = faceOffset + firstInRow(div + 1)
			} else if (faceB.edges()[1].index == edgeIdx) {
				neighbors[1] = faceOffset + lastInRow(subDivs - div)
				neighbors[2] = faceOffset + lastInRow(subDivs - div - 1)
			} else {
				// parallel to rows
				neighbors[1] = faceOffset + lastInRow(subDivs - 2) + div
				neighbors[2] = faceOffset + lastInRow(subDivs - 2) + div + 1
			}
		}

		return neighbors
	}

	vertexNeighborsFace(idx: number, subDivs: number): number[] {
		const [faceIdx, loc] = divmod(idx, (subDivs - 1) * subDivs / 2)
		let neighbors = Array<number>(6)

		// zero indexed row
		const row = Math.ceil((Math.sqrt(8 * (loc + 1) + 1) - 1) * 0.5) - 1

		const along = loc - row * (row + 1) / 2

		const faceOffset = this.baseVerts + this.baseEdges * subDivs + (subDivs - 1) * subDivs / 2 * faceIdx

		// write all in, some will be incorrect
		// row above
		neighbors[0] = faceOffset + loc - row - 1
		neighbors[1] = faceOffset + loc - row
		// next in row
		neighbors[2] = faceOffset + loc + 1
		// row below
		neighbors[3] = faceOffset + loc + row + 2
		neighbors[4] = faceOffset + loc + row + 1
		// prev in row
		neighbors[5] = faceOffset + loc - 1

		// check if first in row, 0 and 5
		if (along == 0) {
			neighbors[0] = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(faceIdx, 0, subDivs - row - 1, subDivs)
			neighbors[5] = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(faceIdx, 0, subDivs - row - 2, subDivs)
		}

		// check if last in row, 1 and 2
		if (along == row) {
			neighbors[1] = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(faceIdx, 1, row, subDivs)
			neighbors[2] = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(faceIdx, 1, row + 1, subDivs)
		}

		// check if in last row, 3 and 4
		if (row == subDivs - 2) {
			neighbors[3] = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(faceIdx, 2, subDivs - along - 2, subDivs)
			neighbors[4] = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(faceIdx, 2, subDivs - along - 1, subDivs)
		}

		return neighbors
	}



	Vertex(idx: number, subDivs: number): Vector {
		const baseVerts = this.baseVerts
		const baseEdges = this.baseEdges
		const baseFaces = this.baseFaces
		// three cases
		// origional vert
		// vert on origional edge
		// vert inside original face
		if (idx < baseVerts) {
			// origional vert
			return this.vertexOrigional(idx, subDivs)
		}
		idx -= baseVerts
		if (idx < baseEdges * subDivs) {
			return this.vertexEdge(idx, subDivs)
		}
		idx -= baseEdges * subDivs
		if (idx < (subDivs - 1) * subDivs / 2 * baseFaces) {
			return this.vertexFace(idx, subDivs)
		}

		throw ("past end")
	}

	vertexOrigional(idx: number, subDivs: number): Vector {
		const coords = this.baseGrid.vertices(idx).coords()
		return coords
	}

	vertexEdge(idx: number, subDivs: number): Vector {
		const [edgeIdx, div] = divmod(idx, subDivs)
		const edge = this.baseGrid.edges(edgeIdx)
		// calcIdx := len(this.baseGrid.Vertices) + edgeIdx*subDivs + div
		// log.Printf("Calculating %d vertex %d divisions along edge %d ", calcIdx, div, edgeIdx)

		const firstVertex = this.baseGrid.vertices(edge.firstVertexA())
		const secondVertex = this.baseGrid.vertices(edge.firstVertexB())
		// angle to subdivide
		let angleToSubdivide = vectorAngle(firstVertex.coords(), secondVertex.coords())

		// angle between origin, first vertex, and second vertex
		let vectorA: [number, number, number] = [0,0,0], vectorB: [number, number, number] = [0,0,0]
		vectorA[0] = -1 * firstVertex.coords()[0]
		vectorA[1] = -1 * firstVertex.coords()[1]
		vectorA[2] = -1 * firstVertex.coords()[2]

		vectorB[0] = secondVertex.coords()[0] - firstVertex.coords()[0]
		vectorB[1] = secondVertex.coords()[1] - firstVertex.coords()[1]
		vectorB[2] = secondVertex.coords()[2] - firstVertex.coords()[2]

		const cornerAngle = vectorAngle(vectorA, vectorB)

		// unit vector from first to second vertex
		let stepDirection: [number, number, number] = [0,0,0]
		stepDirection[0] = vectorB[0] / vectorLength(vectorB)
		stepDirection[1] = vectorB[1] / vectorLength(vectorB)
		stepDirection[2] = vectorB[2] / vectorLength(vectorB)

		// origional radius of the
		const sphereRadius = vectorLength(firstVertex.coords())

		let vert: Vector
		const divisionLength = Math.sin(angleToSubdivide * ((div + 1) / (subDivs + 1))) * sphereRadius / Math.sin(Math.PI - cornerAngle - angleToSubdivide * ((div + 1) / (subDivs + 1)))

		vert[0] = firstVertex.coords()[0] + stepDirection[0] * divisionLength
		vert[1] = firstVertex.coords()[1] + stepDirection[1] * divisionLength
		vert[2] = firstVertex.coords()[2] + stepDirection[2] * divisionLength

		return vert
	}

	vertexFace(idx: number, subDivs: number): Vector {
		const [face, loc] = divmod(idx, (subDivs - 1) * subDivs / 2)

		// log.Printf("Calculating vertex %d locations along division of face %d ", loc, face)

		// zero indexed row
		const row = Math.ceil((Math.sqrt((8 * (loc + 1)) + 1) - 1) * 0.5) - 1

		// calculate trinagle number for previous row [row*(row+1)] since we are zero indexed
		const along = loc - row * (row + 1) / 2

		// get edge verts to work from
		const baseVerts = this.baseVerts
		const first = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(face, 0, subDivs - 2 - row, subDivs - baseVerts)
		const firstVertex = this.vertexEdge(first, subDivs)
		const second = this.baseGrid.vertexIndexAtClockwiseIndexOnOldFace(face, 1, 1 + row, subDivs - baseVerts)
		const secondVertex = this.vertexEdge(second, subDivs)

		// log.Printf("First Vert idx: %d, Second Vert idx: %d", first, second)
		// log.Printf("First %#v, Second: %#v", firstVertex, secondVertex)

		// angle to subdivide
		const angleToSubdivide = vectorAngle(firstVertex, secondVertex)

		// angle between origin, first vertex, and second vertex
		let vectorA: Vector = [0,0,0], vectorB: Vector = [0,0,0]
		vectorA[0] = -1 * firstVertex[0]
		vectorA[1] = -1 * firstVertex[1]
		vectorA[2] = -1 * firstVertex[2]

		vectorB[0] = secondVertex[0] - firstVertex[0]
		vectorB[1] = secondVertex[1] - firstVertex[1]
		vectorB[2] = secondVertex[2] - firstVertex[2]

		const cornerAngle = vectorAngle(vectorA, vectorB)

		// unit vector from first to second vertex
		let stepDirection: Vector = [0,0,0]
		stepDirection[0] = vectorB[0] / vectorLength(vectorB)
		stepDirection[1] = vectorB[1] / vectorLength(vectorB)
		stepDirection[2] = vectorB[2] / vectorLength(vectorB)

		const firstVectorLength = vectorLength(firstVertex)

		let vert: Vector

		const divisionLength = Math.sin(angleToSubdivide * ((along + 1) / (row + 2))) * firstVectorLength / Math.sin(Math.PI - cornerAngle - angleToSubdivide * ((along + 1) / (row + 2)))
		vert[0] = firstVertex[0] + stepDirection[0] * divisionLength
		vert[1] = firstVertex[1] + stepDirection[1] * divisionLength
		vert[2] = firstVertex[2] + stepDirection[2] * divisionLength

		return vert
	}
}

function divmod(numerator: number, denominator: number): [number, number] {
	const quotient = Math.floor(numerator / denominator) // integer division
	const remainder = numerator % denominator
	return [quotient, remainder]
}

function firstInRow(row: number): number {
	return (row - 1) * row / 2
}

function lastInRow(row: number): number {
	return firstInRow(row + 1) - 1
}