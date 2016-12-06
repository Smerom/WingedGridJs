class WingedEdge {
	firstVertexA: number;
	firstVertexB: number;
	faceA: number;
	faceB: number;
	prevA: number;
	nextA: number;
	prevB: number;
	nextB: number;
	/********* Base **********/
	nextEdgeForFace(faceIndex: number): [number, Error] {
		if (this.faceA == faceIndex) {
			return [this.nextA, null];
		} else if (this.faceB == faceIndex){
			return [this.nextB, null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	prevEdgeForFace(faceIndex: number): [number, Error] {
		if (this.faceA == faceIndex) {
			return [this.prevA, null];
		} else if (this.faceB == faceIndex){
			return [this.prevB, null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	firstVertexForFace(faceIndex: number): [number, Error] {
		if (this.faceA == faceIndex) {
			return [this.firstVertexA, null];
		} else if (this.faceB == faceIndex){
			return [this.firstVertexB, null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	secondVertexForFace(faceIndex: number): [number, Error] {
		if (this.faceA == faceIndex) {
			return [this.firstVertexB, null];
		} else if (this.faceB == faceIndex){
			return [this.firstVertexA, null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	nextEdgeForVertex(vertexIndex: number): [number, Error] {
		if (this.firstVertexA == vertexIndex) {
			return this.prevEdgeForFace(this.faceA);
		} else if (this.firstVertexB == vertexIndex){
			return this.prevEdgeForFace(this.faceB);
		}
		return [-1, Error("Edge not associated with vertex.")];
	}
	prevEdgeForVertex(vertexIndex: number): [number, Error] {
		if (this.firstVertexA == vertexIndex) {
			return this.nextEdgeForFace(this.faceB);
		} else if (this.firstVertexB == vertexIndex){
			return this.nextEdgeForFace(this.faceA);
		}
		return [-1, Error("Edge not associated with vertex.")];
	}
	adjacentForFace(faceIndex: number): [number, Error] {
		if (this.faceA == faceIndex) {
			return [this.faceB, null];
		} else if (this.faceB == faceIndex){
			return [this.faceA, null];
		}
		return [-1, Error("Edge not associated with face.")];
	}
	adjacentForVertex(vertexIndex: number): [number, Error] {
		if (this.firstVertexA == vertexIndex) {
			return [this.firstVertexB, null]
		} else if (this.firstVertexB == vertexIndex) {
			return [this.firstVertexA, null]
		}
		return [-1, Error("Edge not associated with vertex.")];
	}
}

class WingedFace {
	edges: Array<number>;
	constructor(){
		this.edges = new Array<number>();
	}
}

class WingedVertex {
	coords: [number, number, number];
	edges: Array<number>;
	vertexNeighbors: Array<number>;
	constructor(){
		this.coords = [0,0,0];
		this.edges = new Array<number>();
		this.vertexNeighbors = new Array<number>();
	}
}

class WingedGrid {
	faces: Array<WingedFace>;
	edges: Array<WingedEdge>;
	vertices: Array<WingedVertex>;
	constructor(faceCount: number, edgeCount: number, vertCount: number){
		this.faces = new Array<WingedFace>(faceCount);
		for (var i = this.faces.length - 1; i >= 0; i--) {
			this.faces[i] = new WingedFace();
		}
		this.edges = new Array<WingedEdge>(edgeCount);
		for (var i = this.edges.length - 1; i >= 0; i--) {
			this.edges[i] = new WingedEdge();
		}
		this.vertices = new Array<WingedVertex>(vertCount);
		for (var i = this.vertices.length - 1; i >= 0; i--) {
			this.vertices[i] = new WingedVertex();
		}
	}
	/********* Base **********/
	neighborsForFace(faceIndex: number): Array<number>{
		let neighbors = Array<number>();
		for (let i = this.faces[faceIndex].edges.length - 1; i >= 0; i--) {
			let [neighbor, err] = this.edges[this.faces[faceIndex].edges[i]].adjacentForFace(faceIndex);
			neighbors[i] = neighbor;
		}
		return neighbors;
	}
	neighborsForVertex(vertexIndex: number): Array<number> {
		let neighbors = Array<number>();
		if(this.vertices[vertexIndex].vertexNeighbors.length == 0) {
			for (var i = this.vertices[vertexIndex].edges.length - 1; i >= 0; i--) {
				let [neighbor, error] = this.edges[this.vertices[vertexIndex].edges[i]].adjacentForVertex(vertexIndex);
				this.vertices[vertexIndex].vertexNeighbors[i] = neighbor;
			}
		}
		return neighbors
	}

	/*********** SPHERE ************/
	normalizeVerticesToDistanceFromOrigin(wantedLength: number) {
		for (var i = this.vertices.length - 1; i >= 0; i--) {
			let vertex = this.vertices[i];
			let currentLength = vectorLength(vertex.coords);
			this.vertices[i].coords[0] = vertex.coords[0] * wantedLength / currentLength;
			this.vertices[i].coords[1] = vertex.coords[1] * wantedLength / currentLength;
			this.vertices[i].coords[2] = vertex.coords[2] * wantedLength / currentLength;
		}
	}

	/******** DUAL *********/
	createDual(): [WingedGrid, Error]{
		let dualGrid = new WingedGrid(this.vertices.length, this.edges.length, this.faces.length); // swap vert and face counts
		// create faces
		for (var i = this.vertices.length - 1; i >= 0; i--) {
			let vertex = this.vertices[i];
			for (var j = vertex.edges.length - 1; j >= 0; j--) {
				dualGrid.faces[i].edges[j] = vertex.edges[j];
			}
		}

		// create vertices
		for (var i = this.faces.length - 1; i >= 0; i--) {
			let face = this.faces[i];
			// set edges
			for (var j = face.edges.length - 1; j >= 0; j--) {
				dualGrid.vertices[i].edges[j] = face.edges[j];
			}
			// set coords from center of old face
			let faceCenter: [number, number, number] = [0,0,0];
			let count: number = 0;
			for (var j = face.edges.length - 1; j >= 0; j--) {
				let edge = this.edges[face.edges[j]];
				let [vertexIndex, err] = edge.firstVertexForFace(i);
				if (err != null) {
					return [dualGrid, err]
				}
				let vertex = this.vertices[vertexIndex];
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
			let edge = this.edges[i];
			dualGrid.edges[i].faceA = edge.firstVertexA;
			dualGrid.edges[i].faceB = edge.firstVertexB;
			// faces are swapped from verts
			dualGrid.edges[i].firstVertexA = edge.faceB;
			dualGrid.edges[i].firstVertexB = edge.faceA;
		}

		// set prev and next
		for (let faceIndex = dualGrid.faces.length - 1; faceIndex >= 0; faceIndex--) {
			let face = dualGrid.faces[faceIndex];
			for (var faceEdgeIndex = face.edges.length - 1; faceEdgeIndex >= 0; faceEdgeIndex--) {
				let edgeIndex = face.edges[faceEdgeIndex];
				let theEdge = dualGrid.edges[edgeIndex];
				if (theEdge.faceA == faceIndex) {
					if (faceEdgeIndex == 0) {
						dualGrid.edges[edgeIndex].prevA = face.edges[face.edges.length - 1];
						dualGrid.edges[edgeIndex].nextA = face.edges[1];
					} else if (faceEdgeIndex == face.edges.length - 1) {
						dualGrid.edges[edgeIndex].prevA = face.edges[face.edges.length - 1];
						dualGrid.edges[edgeIndex].nextA = face.edges[0];
					} else {
						dualGrid.edges[edgeIndex].prevA = face.edges[faceEdgeIndex - 1];
						dualGrid.edges[edgeIndex].nextA = face.edges[faceEdgeIndex + 1];
					}
				} 
				if (theEdge.faceB == faceIndex) {
					if (faceEdgeIndex == 0) {
						dualGrid.edges[edgeIndex].prevB = face.edges[face.edges.length - 1];
						dualGrid.edges[edgeIndex].nextB = face.edges[1];
					} else if (faceEdgeIndex == face.edges.length - 1) {
						dualGrid.edges[edgeIndex].prevB = face.edges[face.edges.length - 1];
						dualGrid.edges[edgeIndex].nextB = face.edges[0];
					} else {
						dualGrid.edges[edgeIndex].prevB = face.edges[faceEdgeIndex - 1];
						dualGrid.edges[edgeIndex].nextB = face.edges[faceEdgeIndex + 1];
					}
				}
			}
		}
		return [dualGrid, null]
	}

	/************* Subdivision helpers ****************/
	// rename vertexInexAtClockwiseIndexOnFaceToSubdivide
	vertexIndexAtClockwiseIndexOnOldFace(faceIndex: number, edgeInFaceIndex: number, clockwiseVertexIndex: number, edgeSubdivisions: number): number {
		let edgeIndex = this.faces[faceIndex].edges[edgeInFaceIndex];
		let edge = this.edges[edgeIndex];
		if (edge.faceA == faceIndex) {
			return this.vertices.length + edgeIndex*edgeSubdivisions + clockwiseVertexIndex;
		}
		if (edge.faceB == faceIndex) {
			return this.vertices.length + edgeIndex*edgeSubdivisions + edgeSubdivisions - 1 - clockwiseVertexIndex;
		}
		return -1;
	}
	edgeIndexAtClockwiseIndexOnOldFace(faceIndex: number, edgeInFaceIndex: number, clockwiseEdgeIndex: number, edgeSubdivisions: number): number {
		let oldEdgeIndex = this.faces[faceIndex].edges[edgeInFaceIndex];
		let oldEdge = this.edges[oldEdgeIndex];
		if (oldEdge.faceA == faceIndex) {
			return oldEdgeIndex*(edgeSubdivisions + 1) + clockwiseEdgeIndex;
		}
		if (oldEdge.faceB == faceIndex) {
			return oldEdgeIndex*(edgeSubdivisions + 1) + edgeSubdivisions - clockwiseEdgeIndex;
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
		const faceCount: number = this.faces.length * ((edgeSubdivisions + 2)*(edgeSubdivisions+1)/2 + (edgeSubdivisions+1)*edgeSubdivisions/2);		
		let dividedGrid = new WingedGrid(
			faceCount,
			// since each face 'owns' 1/2 of three edges, there are 1.5 times as
			//  many edges as faces
			3*faceCount/2,
			// Euler-somebody or other gives us the vertex count of
			//  1/2*(face count) + 2
			faceCount/2 + 2
			);

		// Invalidate all values that will be set later, useful for error checking
		for (let i = 0; i < faceCount; i++) {
			dividedGrid.faces[i].edges = [-1, -1, -1];
		}
		for (let i = 0; i < 3*faceCount/2; i++) {
			dividedGrid.edges[i].faceA = -1
			dividedGrid.edges[i].faceB = -1
			dividedGrid.edges[i].firstVertexA = -1
			dividedGrid.edges[i].firstVertexB = -1
			dividedGrid.edges[i].nextA = -1
			dividedGrid.edges[i].nextB = -1
			dividedGrid.edges[i].prevA = -1
			dividedGrid.edges[i].prevB = -1
		}
		// vertices corrisponding with old ones will have the same number of
		//  associated edges to preserve the Euler characteristic ( =2 for S2)
		for (let i = 0; i < this.vertices.length; i++) {
			dividedGrid.vertices[i].edges = new Array<number>(this.vertices[i].edges.length)
			for (let j = 0; j < this.vertices[i].edges.length; j++) {
				dividedGrid.vertices[i].edges[j] = -1
			}
			// set the coords
			dividedGrid.vertices[i].coords = [Infinity, Infinity, Infinity];
		}
		// the way we divide the faces creates six accociated edges for the remaining
		//  vertecies
		for (let i = this.vertices.length; i < faceCount/2+2; i++) {
			dividedGrid.vertices[i].edges = [-1, -1, -1, -1, -1, -1];
			// set the coords
			dividedGrid.vertices[i].coords = [Infinity, Infinity, Infinity];
		}

		/***************** Subdivide the grid ****************/

		/******************* EDGE SUBDIVISION ***********************/
		/****** Old Edge Subdivision goes in the first section of the new array, ordered by edge *******/
		let origVertexCount = this.vertices.length;
		let origEdgeCount = this.edges.length;
		for (let i = this.edges.length - 1; i >= 0; i--) {
			let edge = this.edges[i]
	
			// first edge has origional vertex
			dividedGrid.edges[i*(edgeSubdivisions+1)].firstVertexA = edge.firstVertexA
			dividedGrid.edges[i*(edgeSubdivisions+1)].firstVertexB = origVertexCount + i*edgeSubdivisions

			for (let j = 1; j < edgeSubdivisions; j++) {
				// set the edge vertex indecies
				dividedGrid.edges[i*(edgeSubdivisions+1)+j].firstVertexA = origVertexCount + i*edgeSubdivisions + j - 1
				dividedGrid.edges[i*(edgeSubdivisions+1)+j].firstVertexB = origVertexCount + i*edgeSubdivisions + j
			}

			// connect last new edge
			dividedGrid.edges[i*(edgeSubdivisions+1)+edgeSubdivisions].firstVertexA = origVertexCount + i*edgeSubdivisions + edgeSubdivisions - 1
			dividedGrid.edges[i*(edgeSubdivisions+1)+edgeSubdivisions].firstVertexB = edge.firstVertexB
		}

		/********* Edges created interior to old faces go in the second section, ordered by face. ****/
		for (let faceIndex = this.faces.length - 1; faceIndex >= 0; faceIndex--) {
			let oldFace = this.faces[faceIndex];
			// vertex offset for vertices interior to the face
			let vertexOffset = origVertexCount + origEdgeCount*edgeSubdivisions + (edgeSubdivisions*(edgeSubdivisions-1)/2)*faceIndex;

			var edgeOffset: number = (edgeSubdivisions+1)*origEdgeCount + 3*(edgeSubdivisions)*(edgeSubdivisions+1)/2*faceIndex;
			if (edgeSubdivisions == 1) {
				dividedGrid.edges[edgeOffset+0].firstVertexA = origVertexCount + oldFace.edges[0];
				dividedGrid.edges[edgeOffset+0].firstVertexB = origVertexCount + oldFace.edges[1];

				dividedGrid.edges[edgeOffset+1].firstVertexA = origVertexCount + oldFace.edges[0];
				dividedGrid.edges[edgeOffset+1].firstVertexB = origVertexCount + oldFace.edges[2];

				dividedGrid.edges[edgeOffset+2].firstVertexA = origVertexCount + oldFace.edges[1];
				dividedGrid.edges[edgeOffset+2].firstVertexB = origVertexCount + oldFace.edges[2];
			} else {
				// first row
				dividedGrid.edges[edgeOffset+0].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions-1, edgeSubdivisions);
				dividedGrid.edges[edgeOffset+0].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 0, edgeSubdivisions);

				dividedGrid.edges[edgeOffset+1].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions-1, edgeSubdivisions);
				dividedGrid.edges[edgeOffset+1].firstVertexB = vertexOffset + 0;

				dividedGrid.edges[edgeOffset+2].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 0, edgeSubdivisions);
				dividedGrid.edges[edgeOffset+2].firstVertexB = vertexOffset + 0;

				// middle rows
				var rowOffset: number;
				var i: number = 0;
				for (i = 1; i < edgeSubdivisions-1; i++) {
					rowOffset = i * (i + 1) * 3 / 2
					// first border
					dividedGrid.edges[edgeOffset+rowOffset+0].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions-1-i, edgeSubdivisions);
					dividedGrid.edges[edgeOffset+rowOffset+0].firstVertexB = vertexOffset + (i * (i - 1) / 2);

					dividedGrid.edges[edgeOffset+rowOffset+1].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions-1-i, edgeSubdivisions);
					dividedGrid.edges[edgeOffset+rowOffset+1].firstVertexB = vertexOffset + (i * (i - 1) / 2) + i;

					dividedGrid.edges[edgeOffset+rowOffset+2].firstVertexA = vertexOffset + (i * (i - 1) / 2);
					dividedGrid.edges[edgeOffset+rowOffset+2].firstVertexB = vertexOffset + (i*(i-1))/2 + i;
					// interior of face
					for (let j = 1; j < i; j++) {
						dividedGrid.edges[edgeOffset+rowOffset+j*3+0].firstVertexA = vertexOffset + (i * (i - 1) / 2) + j - 1;
						dividedGrid.edges[edgeOffset+rowOffset+j*3+0].firstVertexB = vertexOffset + (i * (i - 1) / 2) + j;

						dividedGrid.edges[edgeOffset+rowOffset+j*3+1].firstVertexA = vertexOffset + (i * (i - 1) / 2) + j - 1;
						dividedGrid.edges[edgeOffset+rowOffset+j*3+1].firstVertexB = vertexOffset + (i * (i - 1) / 2) + i + j;

						dividedGrid.edges[edgeOffset+rowOffset+j*3+2].firstVertexA = vertexOffset + (i * (i - 1) / 2) + j;
						dividedGrid.edges[edgeOffset+rowOffset+j*3+2].firstVertexB = vertexOffset + (i * (i - 1) / 2) + i + j;
					}

					// second border
					dividedGrid.edges[edgeOffset+rowOffset+i*3+0].firstVertexA = vertexOffset + (i * (i - 1) / 2) + i - 1;
					dividedGrid.edges[edgeOffset+rowOffset+i*3+0].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);

					dividedGrid.edges[edgeOffset+rowOffset+i*3+1].firstVertexA = vertexOffset + (i * (i - 1) / 2) + i - 1;
					dividedGrid.edges[edgeOffset+rowOffset+i*3+1].firstVertexB = vertexOffset + (i * (i - 1) / 2) + 2*i;

					dividedGrid.edges[edgeOffset+rowOffset+i*3+2].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);
					dividedGrid.edges[edgeOffset+rowOffset+i*3+2].firstVertexB = vertexOffset + (i * (i - 1) / 2) + 2*i;
				}

				// last row
				i = edgeSubdivisions - 1; // should already be set, but just incase
				rowOffset = i * (i + 1) * 3 / 2;
				// border 1
				dividedGrid.edges[edgeOffset+rowOffset+0].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, 0, edgeSubdivisions);
				dividedGrid.edges[edgeOffset+rowOffset+0].firstVertexB = vertexOffset + (i * (i - 1) / 2);

				dividedGrid.edges[edgeOffset+rowOffset+1].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, 0, edgeSubdivisions);
				dividedGrid.edges[edgeOffset+rowOffset+1].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions-1, edgeSubdivisions);

				dividedGrid.edges[edgeOffset+rowOffset+2].firstVertexA = vertexOffset + (i * (i - 1) / 2);
				dividedGrid.edges[edgeOffset+rowOffset+2].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions-1, edgeSubdivisions);

				// middle
				for (let j = 1; j < i; j++) {
					dividedGrid.edges[edgeOffset+rowOffset+j*3+0].firstVertexA = vertexOffset + (i * (i - 1) / 2) + j - 1;
					dividedGrid.edges[edgeOffset+rowOffset+j*3+0].firstVertexB = vertexOffset + (i * (i - 1) / 2) + j;

					dividedGrid.edges[edgeOffset+rowOffset+j*3+1].firstVertexA = vertexOffset + (i * (i - 1) / 2) + j - 1;
					dividedGrid.edges[edgeOffset+rowOffset+j*3+1].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions-j-1, edgeSubdivisions);

					dividedGrid.edges[edgeOffset+rowOffset+j*3+2].firstVertexA = vertexOffset + (i * (i - 1) / 2) + j;
					dividedGrid.edges[edgeOffset+rowOffset+j*3+2].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions-j-1, edgeSubdivisions);
				}

				// border 2
				dividedGrid.edges[edgeOffset+rowOffset+i*3+0].firstVertexA = vertexOffset + (i * (i - 1) / 2) + i - 1;
				dividedGrid.edges[edgeOffset+rowOffset+i*3+0].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);

				dividedGrid.edges[edgeOffset+rowOffset+i*3+1].firstVertexA = vertexOffset + (i * (i - 1) / 2) + i - 1;
				dividedGrid.edges[edgeOffset+rowOffset+i*3+1].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, 0, edgeSubdivisions);

				dividedGrid.edges[edgeOffset+rowOffset+i*3+2].firstVertexA = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions);
				dividedGrid.edges[edgeOffset+rowOffset+i*3+2].firstVertexB = this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 2, 0, edgeSubdivisions);
			}
		}

		// create the faces, update the edges
		/******************* FACE SUBDIVISION ***********************/
		for (let faceIndex = this.faces.length - 1; faceIndex >= 0; faceIndex--) {
			// get the number of faces we divide this one into
			let subFaceCount = (edgeSubdivisions+2)*(edgeSubdivisions+1)/2 + (edgeSubdivisions+1)*edgeSubdivisions/2;

			// number of internal edges created on this face
			//  equal to ( (edgeSubdivisions + 1) choose 2 ) * 3
			var subEdgeCount = 3 * edgeSubdivisions * (edgeSubdivisions + 1) / 2;

			// each new face set starts at
			var indexStart: number = subFaceCount * faceIndex;

			var edgeOffset: number = (edgeSubdivisions+1)*origEdgeCount + subEdgeCount*faceIndex;

			// first corner
			// set edges in face
			dividedGrid.faces[indexStart+0].edges = 
				[ this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions, edgeSubdivisions)
				, this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 0, edgeSubdivisions)
				, edgeOffset];

			// loop through the middle section of faces
			// edges grow by 1/2*i*(i-1)*3
			for (let i = 1; i < edgeSubdivisions; i++) {
				var rowIndexStart: number = i*(i+1)/2 + i*(i-1)/2

				// edge
				dividedGrid.faces[indexStart+rowIndexStart+0].edges = 
					[ this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions-i, edgeSubdivisions)
					, edgeOffset + i*(i-1)*3/2 + 1
					, edgeOffset + i*(i+1)*3/2];

				// middle
				// up to one less than next index start
				for (let j = 1; j < (i+1)*(i+2)/2+(i+1)*i/2-rowIndexStart-1; j++) {
					if (j%2 == 1) {
						dividedGrid.faces[indexStart+rowIndexStart+j].edges = 
							[ edgeOffset + i*(i-1)*3/2 + (j-1)*3/2
							, edgeOffset + i*(i-1)*3/2 + (j-1)*3/2 + 2
							, edgeOffset + i*(i-1)*3/2 + (j-1)*3/2 + 1];
					} else {
						dividedGrid.faces[indexStart+rowIndexStart+j].edges = 
							[ edgeOffset + i*(i+1)*3/2 + j*3/2
							, edgeOffset + i*(i-1)*3/2 + (j-2)*3/2 + 2
							, edgeOffset + i*(i-1)*3/2 + j*3/2 + 1];
					}
				}

				// edge
				dividedGrid.faces[indexStart+(i+1)*(i+2)/2+(i+1)*i/2-1].edges =
					[ edgeOffset + i*(i+1)*3/2 + ((i+1)*(i+2)/2+(i+1)*i/2-rowIndexStart-1)*3/2
					, edgeOffset + i*(i-1)*3/2 + ((i+1)*(i+2)/2+(i+1)*i/2-rowIndexStart-3)*3/2 + 2
					, this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 1, i, edgeSubdivisions)];
			}

			var rowIndexStart = edgeSubdivisions*(edgeSubdivisions+1)/2 + edgeSubdivisions*(edgeSubdivisions-1)/2
			// bottom corner 1
			dividedGrid.faces[indexStart+rowIndexStart+0].edges =
				[ this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 0, 0, edgeSubdivisions)
				, edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + 1
				, this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions, edgeSubdivisions)];
			// bottom edge
			for (let j = 1; j < subFaceCount-rowIndexStart-1; j++) {
				if (j%2 == 1) {
					dividedGrid.faces[indexStart+rowIndexStart+j].edges =
						[ edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + (j-1)*3/2
						, edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + (j-1)*3/2 + 2
						, edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + (j-1)*3/2 + 1];
				} else {
					dividedGrid.faces[indexStart+rowIndexStart+j].edges =
						[ this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 2, edgeSubdivisions-j/2, edgeSubdivisions)
						, edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + (j-2)*3/2 + 2
						, edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + j*3/2 + 1];
				}
			}
			// bottom corner 2
			dividedGrid.faces[indexStart+subFaceCount-1].edges =
				[ this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 2, 0, edgeSubdivisions)
				, edgeOffset + edgeSubdivisions*(edgeSubdivisions-1)*3/2 + (subFaceCount-rowIndexStart-3)*3/2 + 2
				, this.edgeIndexAtClockwiseIndexOnOldFace(faceIndex, 1, edgeSubdivisions, edgeSubdivisions)];
		}

		// set edge faces from the previously build edge arrays
		for (var faceIndex = dividedGrid.faces.length - 1; faceIndex >= 0; faceIndex--) {
			let edges = dividedGrid.faces[faceIndex].edges;
			let edgesLength = edges.length;
			var thisEdge: WingedEdge, nextEdge: WingedEdge, prevEdge: WingedEdge;

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
			let prevEdge = dividedGrid.edges[edges[edgesLength-1]]
			let thisEdge = dividedGrid.edges[edges[0]]
			let nextEdge = dividedGrid.edges[edges[1]]
			// test vertices
			if (thisEdge.firstVertexA == prevEdge.firstVertexA || thisEdge.firstVertexA == prevEdge.firstVertexB) {
				// check the next edge also matches and face A is not set
				if (thisEdge.firstVertexB == nextEdge.firstVertexA || thisEdge.firstVertexB == nextEdge.firstVertexB) {
					if (thisEdge.faceA == -1) {
						dividedGrid.edges[edges[0]].faceA = faceIndex
						dividedGrid.edges[edges[0]].prevA = edges[edgesLength-1]
						dividedGrid.edges[edges[0]].nextA = edges[1]
					} else {
						console.log("For face "+faceIndex+". Face A has already been set for edge: "+thisEdge+" With edge set: "+edges);
					}

				} else {
					console.log("For face "+faceIndex+". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else if (thisEdge.firstVertexB == prevEdge.firstVertexA || thisEdge.firstVertexB == prevEdge.firstVertexB) {
				// check the next edge also matches and face B is not set
				if (thisEdge.firstVertexA == nextEdge.firstVertexA || thisEdge.firstVertexA == nextEdge.firstVertexB) {
					if (thisEdge.faceB == -1) {
						dividedGrid.edges[edges[0]].faceB = faceIndex
						dividedGrid.edges[edges[0]].prevB = edges[edgesLength-1]
						dividedGrid.edges[edges[0]].nextB = edges[1]
					} else {
						console.log("For face "+faceIndex+". Face B has already been set for edge: "+thisEdge+" With edge set: "+edges);
					}
				} else {
					console.log("For face "+faceIndex+". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else {
				console.log("For face "+faceIndex+". Edges Don't share a vertex!");
			}
			
			// loop through middle edges
			for (let i = 1; i < edgesLength-1; i++) {
				prevEdge = dividedGrid.edges[edges[i-1]]
				thisEdge = dividedGrid.edges[edges[i]]
				nextEdge = dividedGrid.edges[edges[i+1]]
				// test vertcies
				if (thisEdge.firstVertexA == prevEdge.firstVertexA || thisEdge.firstVertexA == prevEdge.firstVertexB) {
					// check the next edge also matches
					if (thisEdge.firstVertexB == nextEdge.firstVertexA || thisEdge.firstVertexB == nextEdge.firstVertexB) {
						if (thisEdge.faceA == -1) {
							dividedGrid.edges[edges[i]].faceA = faceIndex
							dividedGrid.edges[edges[i]].prevA = edges[i-1]
							dividedGrid.edges[edges[i]].nextA = edges[i+1]
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
							dividedGrid.edges[edges[i]].faceB = faceIndex
							dividedGrid.edges[edges[i]].prevB = edges[i-1]
							dividedGrid.edges[edges[i]].nextB = edges[i+1]
						} else {
							console.log("For face "+faceIndex+". Face B has already been set for edge: "+thisEdge+" With edge set: "+edges);
						}

					} else {
						console.log("For face "+faceIndex+". Previous edge matches, but next edge doesn't share correct vertex!");
					}
				} else {
					console.log("For face "+faceIndex+". Edges Don't share a vertex!");
				}
			}

			// last edge
			prevEdge = dividedGrid.edges[edges[edgesLength-2]]
			thisEdge = dividedGrid.edges[edges[edgesLength-1]]
			nextEdge = dividedGrid.edges[edges[0]]
			// test vertices
			if (thisEdge.firstVertexA == prevEdge.firstVertexA || thisEdge.firstVertexA == prevEdge.firstVertexB) {
				// check the next edge also matches
				if (thisEdge.firstVertexB == nextEdge.firstVertexA || thisEdge.firstVertexB == nextEdge.firstVertexB) {
					if (thisEdge.faceA == -1) {
						dividedGrid.edges[edges[edgesLength-1]].faceA = faceIndex;
						dividedGrid.edges[edges[edgesLength-1]].prevA = edges[edgesLength-2];
						dividedGrid.edges[edges[edgesLength-1]].nextA = edges[0];
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
						dividedGrid.edges[edges[edgesLength-1]].faceB = faceIndex
						dividedGrid.edges[edges[edgesLength-1]].prevB = edges[edgesLength-2]
						dividedGrid.edges[edges[edgesLength-1]].nextB = edges[0]
					} else {
						console.log("For face "+faceIndex+". Face B has already been set for edge: "+thisEdge+" With edge set: "+edges);
					}
				} else {
					console.log("For face "+faceIndex+". Previous edge matches, but next edge doesn't share correct vertex!");
				}
			} else {
				console.log("For face "+faceIndex+". Edges Don't share a vertex!");
			}//*/
		}

		// create the vertices
		/******************* VERTEX SUBDIVISION ***********************/
		// set coords for the origional verts
		for (let i = this.vertices.length - 1; i >= 0; i--) {
			const vertex = this.vertices[i];
			dividedGrid.vertices[i].coords[0] = vertex.coords[0]
			dividedGrid.vertices[i].coords[1] = vertex.coords[1]
			dividedGrid.vertices[i].coords[2] = vertex.coords[2]
		}

		// subdivide along each edge
		for (let i = this.edges.length - 1; i >= 0; i--) {
			let edge = this.edges[i];
			let firstVertex = this.vertices[edge.firstVertexA];
			let secondVertex = this.vertices[edge.firstVertexB];
			// angle to subdivide
			let angleToSubdivide: number;
			angleToSubdivide = vectorAngle(firstVertex.coords, secondVertex.coords)

			// angle between origin, first vertex, and second vertex
			let vectorA: [number, number, number] = [-1, -1, -1], vectorB: [number, number, number] = [-1,-1,-1];
			vectorA[0] = -1 * firstVertex.coords[0];
			vectorA[1] = -1 * firstVertex.coords[1];
			vectorA[2] = -1 * firstVertex.coords[2];

			vectorB[0] = secondVertex.coords[0] - firstVertex.coords[0];
			vectorB[1] = secondVertex.coords[1] - firstVertex.coords[1];
			vectorB[2] = secondVertex.coords[2] - firstVertex.coords[2];

			let cornerAngle: number;
			cornerAngle = vectorAngle(vectorA, vectorB);

			// unit vector from first to second vertex
			let stepDirection: [number, number, number] = [-1,-1,-1];
			stepDirection[0] = vectorB[0] / vectorLength(vectorB);
			stepDirection[1] = vectorB[1] / vectorLength(vectorB);
			stepDirection[2] = vectorB[2] / vectorLength(vectorB);

			// origional radius of the
			let sphereRadius: number = vectorLength(firstVertex.coords);

			let divisionLength: number;
			for (let j = 0; j < edgeSubdivisions; j++) {
				// find the new vertex position and create the vertex
				// but don't correct it's length yet
				divisionLength = Math.sin(angleToSubdivide*((j+1)/(edgeSubdivisions+1))) * sphereRadius / Math.sin(Math.PI-cornerAngle-angleToSubdivide*((j+1)/(edgeSubdivisions+1)));

				dividedGrid.vertices[origVertexCount+i*edgeSubdivisions+j].coords[0] = firstVertex.coords[0] + stepDirection[0]*divisionLength;
				dividedGrid.vertices[origVertexCount+i*edgeSubdivisions+j].coords[1] = firstVertex.coords[1] + stepDirection[1]*divisionLength;
				dividedGrid.vertices[origVertexCount+i*edgeSubdivisions+j].coords[2] = firstVertex.coords[2] + stepDirection[2]*divisionLength;
			}
		}

		// subdivide face interior
		for (let faceIndex = 0; faceIndex < this.faces.length; faceIndex++) {
			var vertexOffset: number = origVertexCount + origEdgeCount*edgeSubdivisions + edgeSubdivisions*(edgeSubdivisions-1)/2*faceIndex;
			// only if we have more than one division
			if (edgeSubdivisions > 1) {
				for (let i = 0; i < edgeSubdivisions-1; i++) {
					var firstVertex = dividedGrid.vertices[this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 0, edgeSubdivisions-2-i, edgeSubdivisions)];
					var secondVertex = dividedGrid.vertices[this.vertexIndexAtClockwiseIndexOnOldFace(faceIndex, 1, 1+i, edgeSubdivisions)];
					// angle to subdivide
					var angleToSubdivide: number;
					angleToSubdivide = vectorAngle(firstVertex.coords, secondVertex.coords);

					// angle between origin, first vertex, and second vertex
					var vectorA: [number, number, number] = [-1,-1,-1], vectorB: [number, number, number] = [-1,-1,-1];
					vectorA[0] = -1 * firstVertex.coords[0];
					vectorA[1] = -1 * firstVertex.coords[1];
					vectorA[2] = -1 * firstVertex.coords[2];

					vectorB[0] = secondVertex.coords[0] - firstVertex.coords[0];
					vectorB[1] = secondVertex.coords[1] - firstVertex.coords[1];
					vectorB[2] = secondVertex.coords[2] - firstVertex.coords[2];

					var cornerAngle: number;
					cornerAngle = vectorAngle(vectorA, vectorB);

					// unit vector from first to second vertex
					var stepDirection: [number, number, number] = [-1,-1,-1];
					stepDirection[0] = vectorB[0] / vectorLength(vectorB);
					stepDirection[1] = vectorB[1] / vectorLength(vectorB);
					stepDirection[2] = vectorB[2] / vectorLength(vectorB);

					var firstVectorLength: number = vectorLength(firstVertex.coords);
					var divisionLength: number = 0
					for (let j = 0; j < i+1; j++) {
						divisionLength = Math.sin(angleToSubdivide*((j+1)/(i+2))) * firstVectorLength / Math.sin(Math.PI-cornerAngle-angleToSubdivide*((j+1)/(i+2)));
						if (vertexOffset+(i*(i+1)/2)+j >= dividedGrid.vertices.length) {
							console.log("breakpoint");
						}
						dividedGrid.vertices[vertexOffset+(i*(i+1)/2)+j].coords[0] = firstVertex.coords[0] + stepDirection[0]*divisionLength;
						dividedGrid.vertices[vertexOffset+(i*(i+1)/2)+j].coords[1] = firstVertex.coords[1] + stepDirection[1]*divisionLength;
						dividedGrid.vertices[vertexOffset+(i*(i+1)/2)+j].coords[2] = firstVertex.coords[2] + stepDirection[2]*divisionLength;
					}
				}
			}
		}

		// set vertex edge array
		// loop through edges so we only have to touch each one once
		for (let index = dividedGrid.edges.length - 1; index >= 0; index--) {
			var edge = dividedGrid.edges[index];
			var nextEdgeIndex: number = -1;
			var nextEdge: WingedEdge;
			var theVertexIndex: number = edge.firstVertexA;
			var theVertex: WingedVertex = dividedGrid.vertices[theVertexIndex];
			var error: Error;
			if (theVertex.edges[0] == -1) {
				[nextEdgeIndex, error] = edge.nextEdgeForVertex(theVertexIndex);
				nextEdge = dividedGrid.edges[nextEdgeIndex];
				theVertex.edges[0] = index;
				var i: number = 1;
				while (index != nextEdgeIndex) {
					theVertex.edges[i] = nextEdgeIndex;
					[nextEdgeIndex, error] = nextEdge.nextEdgeForVertex(theVertexIndex);
					nextEdge = dividedGrid.edges[nextEdgeIndex];
					i = i + 1;
				}
			}
			theVertexIndex = edge.firstVertexB;
			theVertex = dividedGrid.vertices[theVertexIndex];
			if (theVertex.edges[0] == -1) {
				[nextEdgeIndex, error] = edge.nextEdgeForVertex(theVertexIndex);
				nextEdge = dividedGrid.edges[nextEdgeIndex];
				theVertex.edges[0] = index
				var i: number = 1;
				while (index != nextEdgeIndex) {
					theVertex.edges[i] = nextEdgeIndex;
					[nextEdgeIndex, error] = nextEdge.nextEdgeForVertex(theVertexIndex);
					nextEdge = dividedGrid.edges[nextEdgeIndex];
					i = i + 1;
				}
			}
		}

		return [dividedGrid, null];
	}
}

function baseIcosahedron(): [WingedGrid, Error] {
	let icosahedron = new WingedGrid(0,0,0); // arrays will be set manually

	icosahedron.faces = [ 
		// cap 1
		  {edges: [9, 4, 0]} // 0
		, {edges: [1, 5, 0]} // 1
		, {edges: [2, 6, 1]} // 2
		, {edges: [3, 7, 2]} // 3
		, {edges: [4, 8, 3]} // 4
		// cap 2
		, {edges: [14, 19, 10]} // 5
		, {edges: [15, 11, 10]} // 6
		, {edges: [16, 12, 11]} // 7
		, {edges: [17, 13, 12]} // 8
		, {edges: [18, 14, 13]} // 9
		// Ring of 10 between caps, clockwise around cap 1
		//   starting with the first face adjacent to face 0
		, {edges: [29, 28, 9]}  // 10
		, {edges: [29, 20, 17]} // 11
		, {edges: [21, 20, 5]}  // 12
		, {edges: [21, 22, 18]} // 13
		, {edges: [23, 22, 6]}  // 14
		, {edges: [23, 24, 19]} // 15
		, {edges: [25, 24, 7]}  // 16
		, {edges: [25, 26, 15]} // 17
		, {edges: [27, 26, 8]}  // 18
		, {edges: [27, 28, 16]} // 19
	];
	const goldenRatio = 1.61803398875
	icosahedron.vertices = [
		// Edges are duplicate info, could be done pragmatically after constucting
		// edges
		// Must be in counter-clockwise order
		// y-z plane rectangle
		  {coords: [0, 1, goldenRatio],
		    edges: [4, 3, 2, 1, 0],
			vertexNeighbors:[]} // 0
		, {coords: [0, 1, -goldenRatio],
		    edges: [19, 24, 25, 15, 10],
			vertexNeighbors:[]} // 1
		, {coords: [0, -1, goldenRatio],
		    edges: [5, 20, 29, 9, 0],
			vertexNeighbors:[]} // 2
		, {coords: [0, -1, -goldenRatio],
		    edges: [11, 12, 13, 14, 10],
			vertexNeighbors:[]} // 3
		// x-z plane rectangle
		, {coords: [-goldenRatio, 0, 1],
		    edges: [6, 22, 21, 5, 1],
			vertexNeighbors:[]} // 4
		, {coords: [goldenRatio, 0, 1],
		    edges: [9, 28, 27, 8, 4],
			vertexNeighbors:[]} // 5
		, {coords: [-goldenRatio, 0, -1],
		    edges: [18, 22, 23, 19, 14],
			vertexNeighbors:[]} // 6
		, {coords: [goldenRatio, 0, -1],
		    edges: [15, 26, 27, 16, 11],
			vertexNeighbors:[]} // 7
		// x-y plane rectangle
		, {coords: [-1, goldenRatio, 0],
		    edges: [7, 24, 23, 6, 2],
			vertexNeighbors:[]} // 8
		, {coords: [-1, -goldenRatio, 0],
		    edges: [17, 20, 21, 18, 13],
			vertexNeighbors:[]} // 9
		, {coords: [1, goldenRatio, 0],
		    edges: [8, 26, 25, 7, 3],
			vertexNeighbors:[]} // 10
		, {coords: [1, -goldenRatio, 0],
		    edges: [16, 28, 29, 17, 12],
			vertexNeighbors:[]} // 11
	];
	icosahedron.edges = [
		// cap 1 around vertex 0
		// 5 spokes starting between face 0 and 1,
		// clockwise around starting with short edge
		// of the rectangle (y-z)
		<WingedEdge> {
			firstVertexA: 0, firstVertexB: 2,
			faceA: 0, faceB: 1,
			prevA: 4, nextA: 9,
			prevB: 5, nextB: 1,
		}, // 0
		<WingedEdge> {
			firstVertexA: 0, firstVertexB: 4,
			faceA: 1, faceB: 2,
			prevA: 0, nextA: 5,
			prevB: 6, nextB: 2,
		}, // 1
		<WingedEdge> {
			firstVertexA: 0, firstVertexB: 8,
			faceA: 2, faceB: 3,
			prevA: 1, nextA: 6,
			prevB: 7, nextB: 3,
		}, // 2
		<WingedEdge> {
			firstVertexA: 0, firstVertexB: 10,
			faceA: 3, faceB: 4,
			prevA: 2, nextA: 7,
			prevB: 8, nextB: 4,
		}, // 3
		<WingedEdge> {
			firstVertexA: 0, firstVertexB: 5,
			faceA: 4, faceB: 0,
			prevA: 3, nextA: 8,
			prevB: 9, nextB: 0,
		}, // 4
		// ring of 5 around the base of the cap
		<WingedEdge> {
			firstVertexA: 4, firstVertexB: 2,
			faceA: 1, faceB: 12,
			prevA: 1, nextA: 0,
			prevB: 20, nextB: 21,
		}, // 5
		<WingedEdge> {
			firstVertexA: 8, firstVertexB: 4,
			faceA: 2, faceB: 14,
			prevA: 2, nextA: 1,
			prevB: 22, nextB: 23,
		}, // 6
		<WingedEdge> {
			firstVertexA: 10, firstVertexB: 8,
			faceA: 3, faceB: 16,
			prevA: 3, nextA: 2,
			prevB: 24, nextB: 25,
		}, // 7
		<WingedEdge> {
			firstVertexA: 5, firstVertexB: 10,
			faceA: 4, faceB: 18,
			prevA: 4, nextA: 3,
			prevB: 26, nextB: 27,
		}, // 8
		<WingedEdge> {
			firstVertexA: 2, firstVertexB: 5,
			faceA: 0, faceB: 10,
			prevA: 0, nextA: 4,
			prevB: 28, nextB: 29,
		}, // 9
		// cap 2 around vertex 3
		// 5 spokes starting between face 5 and 6,
		// counter-clockwise from short edge (y-z) rectangle
		<WingedEdge> {
			firstVertexA: 1, firstVertexB: 3,
			faceA: 5, faceB: 6,
			prevA: 19, nextA: 14,
			prevB: 11, nextB: 15,
		}, // 10
		<WingedEdge> {
			firstVertexA: 7, firstVertexB: 3,
			faceA: 6, faceB: 7,
			prevA: 15, nextA: 10,
			prevB: 12, nextB: 16,
		}, // 11
		<WingedEdge> {
			firstVertexA: 11, firstVertexB: 3,
			faceA: 7, faceB: 8,
			prevA: 16, nextA: 11,
			prevB: 13, nextB: 17,
		}, // 12
		<WingedEdge> {
			firstVertexA: 9, firstVertexB: 3,
			faceA: 8, faceB: 9,
			prevA: 17, nextA: 12,
			prevB: 14, nextB: 18,
		}, // 13
		<WingedEdge> {
			firstVertexA: 6, firstVertexB: 3,
			faceA: 9, faceB: 5,
			prevA: 18, nextA: 13,
			prevB: 10, nextB: 19,
		}, // 14
		// ring of 5 around the base of cap 2
		<WingedEdge> {
			firstVertexA: 1, firstVertexB: 7,
			faceA: 6, faceB: 17,
			prevA: 10, nextA: 11,
			prevB: 26, nextB: 25,
		}, // 15
		<WingedEdge> {
			firstVertexA: 7, firstVertexB: 11,
			faceA: 7, faceB: 19,
			prevA: 11, nextA: 12,
			prevB: 28, nextB: 27,
		}, // 16
		<WingedEdge> {
			firstVertexA: 11, firstVertexB: 9,
			faceA: 8, faceB: 11,
			prevA: 12, nextA: 13,
			prevB: 20, nextB: 29,
		}, // 17
		<WingedEdge> {
			firstVertexA: 9, firstVertexB: 6,
			faceA: 9, faceB: 13,
			prevA: 13, nextA: 14,
			prevB: 22, nextB: 21,
		}, // 18
		<WingedEdge> {
			firstVertexA: 6, firstVertexB: 1,
			faceA: 5, faceB: 15,
			prevA: 14, nextA: 10,
			prevB: 24, nextB: 23,
		}, // 19
		// zig-zag down the middle
		// 10 triangles, 10 new edges
		// starting clockwise from end of edge 0
		<WingedEdge> {
			firstVertexA: 2, firstVertexB: 9,
			faceA: 11, faceB: 12,
			prevA: 29, nextA: 17,
			prevB: 21, nextB: 5,
		}, // 20
		<WingedEdge> {
			firstVertexA: 4, firstVertexB: 9,
			faceA: 12, faceB: 13,
			prevA: 5, nextA: 20,
			prevB: 18, nextB: 22,
		}, // 21
		<WingedEdge> {
			firstVertexA: 4, firstVertexB: 6,
			faceA: 13, faceB: 14,
			prevA: 21, nextA: 18,
			prevB: 23, nextB: 6,
		}, // 22
		<WingedEdge> {
			firstVertexA: 8, firstVertexB: 6,
			faceA: 14, faceB: 15,
			prevA: 6, nextA: 22,
			prevB: 19, nextB: 24,
		}, // 23
		<WingedEdge> {
			firstVertexA: 8, firstVertexB: 1,
			faceA: 15, faceB: 16,
			prevA: 23, nextA: 19,
			prevB: 25, nextB: 7,
		}, // 24
		<WingedEdge> {
			firstVertexA: 10, firstVertexB: 1,
			faceA: 16, faceB: 17,
			prevA: 7, nextA: 24,
			prevB: 15, nextB: 26,
		}, // 25
		<WingedEdge> {
			firstVertexA: 10, firstVertexB: 7,
			faceA: 17, faceB: 18,
			prevA: 25, nextA: 15,
			prevB: 27, nextB: 8,
		}, // 26
		<WingedEdge> {
			firstVertexA: 5, firstVertexB: 7,
			faceA: 18, faceB: 19,
			prevA: 8, nextA: 26,
			prevB: 16, nextB: 28,
		}, // 27
		<WingedEdge> {
			firstVertexA: 5, firstVertexB: 11,
			faceA: 19, faceB: 10,
			prevA: 27, nextA: 16,
			prevB: 29, nextB: 9,
		}, // 28
		<WingedEdge> {
			firstVertexA: 2, firstVertexB: 11,
			faceA: 10, faceB: 11,
			prevA: 9, nextA: 28,
			prevB: 17, nextB: 20,
		} // 29
	]

	return [icosahedron, null]
}

function vectorAngle(first: [number, number, number], second: [number, number, number]): number {
	return Math.acos((first[0]*second[0] + first[1]*second[1] + first[2]*second[2]) / (Math.sqrt(first[0]*first[0]+first[1]*first[1]+first[2]*first[2]) * Math.sqrt(second[0]*second[0]+second[1]*second[1]+second[2]*second[2])));
}

function vectorLength(vector: [number, number, number]): number {
	return Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}