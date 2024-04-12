import * as CG from './transforms.js';
import { Matrix, Vector } from "./matrix.js";

const LEFT =   32; // binary 100000
const RIGHT =  16; // binary 010000
const BOTTOM = 8;  // binary 001000
const TOP =    4;  // binary 000100
const FAR =    2;  // binary 000010
const NEAR =   1;  // binary 000001
const FLOAT_EPSILON = 0.000001;

class Renderer {
    // canvas:              object ({id: __, width: __, height: __})
    // scene:               object (...see description on Canvas)
    constructor(canvas, scene) {
        this.canvas = document.getElementById(canvas.id);
        this.canvas.width = canvas.width;
        this.canvas.height = canvas.height;
        this.ctx = this.canvas.getContext('2d');
        this.scene = this.processScene(scene);
        this.enable_animation = false;  // <-- disabled for easier debugging; enable for animation
        this.start_time = null;
        this.prev_time = null;
    }

    //
    updateTransforms(time, delta_time) {
        // TODO: update any transformations needed for animation
    }

    generateCircle(num_sides, model) {
        let vertices = [];
        for (let i = 0; i < num_sides; i++) {
            if ( top_bot == 1 && i == 0 ) { // generating vertices of bottom face
                h = -h;
            }
            let currentTheta = (i / num_sides) * (2 * Math.PI);
            vertices.push( CG.Vector4( model.center[0] + model.radius * Math.cos(currentTheta),
                                    model.center[1] + h,
                                    model.center[2] + model.radius * Math.sin(currentTheta), 
                                    1 ) );
        }
        console.log('here');
        return vertices;
    }


    generateModel(model, new_model) {
        let center = model.center;
        let vertices = [];
        let edges = [];
        if (model.type == "cube") {
            // get required dimensions
            let w = model.width/2;
            let h = model.height/2;
            let d = model.depth/2;

            // calculate front face
            vertices.push( CG.Vector4(center[0]-w, center[1]+h, center[2]+d, 1) ); // top left
            vertices.push( CG.Vector4(center[0]+w, center[1]+h, center[2]+d, 1) ); // top right
            vertices.push( CG.Vector4(center[0]+w, center[1]-h, center[2]+d, 1) ); // bottom right
            vertices.push( CG.Vector4(center[0]-w, center[1]-h, center[2]+d, 1) ); // bottom left
            edges.push( [0, 1, 2, 3, 0] );

            // calculate back face
            vertices.push( CG.Vector4(center[0]-w, center[1]+h, center[2]-d, 1) ); // top left
            vertices.push( CG.Vector4(center[0]+w, center[1]+h, center[2]-d, 1) ); // top right
            vertices.push( CG.Vector4(center[0]+w, center[1]-h, center[2]-d, 1) ); // bottom right
            vertices.push( CG.Vector4(center[0]-w, center[1]-h, center[2]-d, 1) ); // bottom left
            edges.push( [4, 5, 6, 7, 4] );

            // add in additional edges
            edges.push( [0, 4] );
            edges.push( [1, 5] );
            edges.push( [2, 6] );
            edges.push( [3, 7] );

        } else if (model.type == "cylinder") {
            let h = model.height/2
            let r = model.radius
            
            // generate top and bottom circular faces
            for (let top_bot = 0; top_bot < 2; top_bot++) {
                for (let i = 0; i < model.sides; i++) {
                    if ( top_bot == 1 && i == 0 ) { // generating vertices of bottom face
                        h = -h;
                    }
                    let currentTheta = (i / model.sides) * (2 * Math.PI);
                    vertices.push( CG.Vector4( model.center[0] + model.radius * Math.cos(currentTheta),
                                               model.center[1] + h,
                                               model.center[2] + model.radius * Math.sin(currentTheta), 
                                               1 ) );
                }
            }
            // loop through model.sides to generate these
            edges.push( [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0] );
            edges.push( [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 12] );
            // edges.push( [0,12]);
            // add the edges between the vertices

        } else if (model.type == "sphere") {
            // let slices = model.slices;
            // let stacks = model.stacks;
            // let radius = model.radius;
            // let bottom = model.center[1]-radius;

            // for (let i = 0; i < 4; i++) {
            //     let currentTheta = (i / stacks) * (2 * Math.PI);
            //     vertices.push( CG.Vector4( model.center[0] + model.radius * Math.cos(currentTheta),
            //                                bottom + (model.radius*2)/(stacks+1)*(i+1),
            //                                model.center[2] + model.radius * Math.sin(currentTheta), 
            //                                1 ) );
            // }

            // edges.push( [0, 1]);


        } else if (model.type == "cone") {


        }
        new_model.vertices = vertices;
        new_model.edges = edges;
    }
    

    //
    rotateLeft() {
        let omega = Math.PI/16;
        // translate SRP to origin
        let srp = this.scene.view.srp;

        let translation_matrix = new Matrix(4,4);
        CG.mat4x4Translate(translation_matrix, srp.x*-1, srp.y*-1, srp.z*-1);
        let origin_srp = Matrix.multiply([translation_matrix, CG.Vector4(srp.x, srp.y, srp.z, 1)])
        // console.log("origin:");
        // console.log(origin_srp);
        
        //rotate SRP at origin
        let rotation_matrix = new Matrix(4,4);
        CG.mat4x4RotateY(rotation_matrix, omega);
        console.log(rotation_matrix);
        let rotated_around_origin = Matrix.multiply([rotation_matrix, origin_srp]); // [origin_srp.x, origin_srp.y, origin_srp.z, 1]])
        // console.log(rotated_around_origin);

        //trasnlate SRP back to original location
        let translation_matrix_2 = CG.Vector4(4,1);
        CG.mat4x4Translate(translation_matrix_2, origin_srp[0]*-1, origin_srp[1]*-1, origin_srp[2]*-1)

        this.scene.view.srp = translation_matrix_2;
        this.draw();
    }
    
    //
    rotateRight() {

    }
    
    //
    moveLeft() {
        this.scene.view.prp.x = Number(this.scene.view.prp.x) - 1;
        this.scene.view.srp.x = Number(this.scene.view.srp.x) - 1;
        this.draw();
    }
    
    //
    moveRight() {
        this.scene.view.prp.x = Number(this.scene.view.prp.x) + 1;
        this.scene.view.srp.x = Number(this.scene.view.srp.x) + 1;
        this.draw();
    }
    
    //
    moveBackward() {
        this.scene.view.prp.z = Number(this.scene.view.prp.z) + 1;
        this.scene.view.srp.z = Number(this.scene.view.srp.z) + 1;
        this.draw();
    }
    
    //
    moveForward() {
        this.scene.view.prp.z = Number(this.scene.view.prp.z) - 1;
        this.scene.view.srp.z = Number(this.scene.view.srp.z) - 1;
        this.draw();
    }

    //
    draw() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        //define this matrix outside of the for loops as this doesn't change (yet, it will when we start doing arrow keys, etc?)
        let perspective_matrix = CG.mat4x4Perspective(this.scene.view.prp, this.scene.view.srp, this.scene.view.vup, this.scene.view.clip);
        let mper_matrix = CG.mat4x4MPer();
        let viewport = CG.mat4x4Viewport(this.canvas.width, this.canvas.height);

        // loop through each model 
        for (let i=0; i< this.scene.models.length; i++) {
            //loop through each vertex in the model
            let vertices = [];
            for (let j=0; j<this.scene.models[i].vertices.length; j++){
                let new_vertex = Matrix.multiply([perspective_matrix, this.scene.models[i].vertices[j]]);
                vertices.push(new_vertex);       
            }

            //   * For (each line segment in) each edge
            for (let k=0; k<this.scene.models[i].edges.length; k++) {

                // loop through the vertices in each edge which have number that correspond to the vertices...
                // For each line segment (in each edge)
                let actual_vertices = []; //actual vertices for the edge
                for (let v=0; v<this.scene.models[i].edges[k].length; v++) {
                    // use the vertex numbers in the edge to create a list with the actual vertices
                    let vertex_number = this.scene.models[i].edges[k][v]; //in model i, in the edge k, get the vertex number v
                    let vertex_from_number = vertices[vertex_number]; // using this number, get its corresponding Vector
                    actual_vertices.push(vertex_from_number); // push the actual vertex, a Vector, onto the list
                }

                // the actual_vertices correspond to the edge --> line
                console.log("actual vertices");
                console.log(actual_vertices);

                let clipped_vertices = [];

                //loop through the actual_vertices and clip each line
                for (let p=0; p<actual_vertices.length-1; p++) {
                    let line = {pt0: actual_vertices[i], pt1: actual_vertices[i+1]}

                    let new_line = this.clipLinePerspective(line, this.scene.view.clip[4]);
                    if (new_line != null) {
                        clipped_vertices.push(new_line[0]);
                        clipped_vertices.push(newline[1]);
                    }
                }   




               // line:         object {pt0: Vector4, pt1: Vector4}
                    // z_min:        float (near clipping plane in canonical view volume)
                // Clip each edge 
                // we have to send the vertices corresponding to the edge
                    // Do this part last/later



                // now our vertices_of_edge are actual Vectors
                //console.log(actual_vertices);

                // Project to 2D by multiplying by MPer
                let vertices_of_edge = []
                for (let p=0; p<actual_vertices.length; p++) {
                    let mper_vertex = Matrix.multiply([mper_matrix, actual_vertices[p]]);
                    vertices_of_edge.push(mper_vertex);
                }
                //console.log(vertices_of_edge);

                // Translate to viewport
                let scaled_vertices = [];
                // For each edge, multiply the vertices (Vectors) by the viewport to scale them appropriately
                for (let vtx=0; vtx<vertices_of_edge.length; vtx++) {
                    let scaled_vertex = Matrix.multiply([viewport, vertices_of_edge[vtx]]);
                    scaled_vertices.push(scaled_vertex);
                }

                //now all of our vertices are Vectors that are scaled appropraitely and we can draw lines between vertices
                // draw the line(s)
                // loop through scaled vertices and draw appropriate lines
                for (let svtx=0; svtx<scaled_vertices.length-1; svtx++) {
                    this.drawLine(scaled_vertices[svtx].x/scaled_vertices[svtx].w, scaled_vertices[svtx].y/scaled_vertices[svtx].w, scaled_vertices[svtx+1].x/scaled_vertices[svtx+1].w, scaled_vertices[svtx+1].y/scaled_vertices[svtx+1].w);
                    // draw in cartesian by dividing by w
                }
            }
        }

        // For each model
        //   * For each vertex
        //     * transform endpoints to canonical view volume
        //   * For each line segment in each edge
        //     * clip in 3D
        //     * project to 2D
        //     * translate/scale to viewport (i.e. window)
        //     * draw line
    }

    // Get outcode for a vertex
    // vertex:       Vector4 (transformed vertex in homogeneous coordinates)
    // z_min:        float (near clipping plane in canonical view volume)
    outcodePerspective(vertex, z_min) {
        let outcode = 0;
        if (vertex.x < (vertex.z - FLOAT_EPSILON)) {
            outcode += LEFT;
        }
        else if (vertex.x > (-vertex.z + FLOAT_EPSILON)) {
            outcode += RIGHT;
        }
        if (vertex.y < (vertex.z - FLOAT_EPSILON)) {
            outcode += BOTTOM;
        }
        else if (vertex.y > (-vertex.z + FLOAT_EPSILON)) {
            outcode += TOP;
        }
        if (vertex.z < (-1.0 - FLOAT_EPSILON)) {
            outcode += FAR;
        }
        else if (vertex.z > (z_min + FLOAT_EPSILON)) {
            outcode += NEAR;
        }
        return outcode;
    }

    calcXYZ(t){
        let x = pt0.x + t*change_x;
        let y = pt0.y + t*change_y;
        let z = pt0.z + t*change_z;
        let intersect_pt = {x: x, y:y, z:z};
        return intersect_pt;    }

    calctLeft(pt0, pt1) {
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = ((-1*pt0.x) + pt0.z)/(change_x-change_z);
        this.calcXYZ(t);
    }

    // calculateIntersectionX(pt0, pt1, x_val) {
    //     let intersect_pt = {x: x_val, y: 0};
    
    //     let t = (x_val - pt0.x)/(pt1.x - pt0.x);
    //     y = (1-t)*pt0.y + t*pt1.y;
    //     intersect_pt.y = y;
    
    //     return intersect_pt;
    // }
                
    calculateIntersectionY(pt0, pt1, y_val) {
        let intersect_pt = {x: 0, y: y_val};
    
        let t = (y_val - pt0.y)/(pt1.y - pt0.y);
        x = (1-t)*pt0.x + t*pt1.x;
        intersect_pt.x = x;
    
        return intersect_pt;
    }

    intersectNear(pt0, pt1, z_val) {
        let intersect_pt = {x: 0, y: y_val};
    
        let t = (y_val - pt0.y)/(pt1.y - pt0.y);
        x = (1-t)*pt0.x + t*pt1.x;
        intersect_pt.x = x;
    
        return intersect_pt;
    }



    // Clip line - should either return a new line (with two endpoints inside view volume)
    //             or null (if line is completely outside view volume)
    // line:         object {pt0: Vector4, pt1: Vector4}
    // z_min:        float (near clipping plane in canonical view volume)
    clipLinePerspective(line, z_min) {
        let result = null;
        let p0 = CG.Vector3(line.pt0.x, line.pt0.y, line.pt0.z); 
        let p1 = CG.Vector3(line.pt1.x, line.pt1.y, line.pt1.z);
        let out0 = this.outcodePerspective(p0, z_min);
        let out1 = this.outcodePerspective(p1, z_min);

        let accept = false;
        while (accept == false) {
            if ((out0 | out1) == 0) {
                accept = true;
                result = line;
            } else if ((out0 & out1) != 0) {
                accept = true;
                result = null;
            } else {
                
            }

        }
        
        // TODO: implement clipping here!
        
        return result;
    }

    //
    animate(timestamp) {
        // Get time and delta time for animation
        if (this.start_time === null) {
            this.start_time = timestamp;
            this.prev_time = timestamp;
        }
        let time = timestamp - this.start_time;
        let delta_time = timestamp - this.prev_time;

        // Update transforms for animation
        this.updateTransforms(time, delta_time);

        // Draw slide
        this.draw();

        // Invoke call for next frame in animation
        if (this.enable_animation) {
            window.requestAnimationFrame((ts) => {
                this.animate(ts);
            });
        }

        // Update previous time to current one for next calculation of delta time
        this.prev_time = timestamp;
    }

    //
    updateScene(scene) {
        this.scene = this.processScene(scene);
        if (!this.enable_animation) {
            this.draw();
        }
    }

    //
    processScene(scene) {
        let processed = {
            view: {
                prp: CG.Vector3(scene.view.prp[0], scene.view.prp[1], scene.view.prp[2]),
                srp: CG.Vector3(scene.view.srp[0], scene.view.srp[1], scene.view.srp[2]),
                vup: CG.Vector3(scene.view.vup[0], scene.view.vup[1], scene.view.vup[2]),
                clip: [...scene.view.clip]
            },
            models: []
        };

        for (let i = 0; i < scene.models.length; i++) {
            let model = { type: scene.models[i].type };
            if (model.type === 'generic') {
                model.vertices = [];
                model.edges = JSON.parse(JSON.stringify(scene.models[i].edges));
                for (let j = 0; j < scene.models[i].vertices.length; j++) {
                    model.vertices.push(CG.Vector4(scene.models[i].vertices[j][0],
                                                   scene.models[i].vertices[j][1],
                                                   scene.models[i].vertices[j][2],
                                                   1));
                    if (scene.models[i].hasOwnProperty('animation')) {
                        model.animation = JSON.parse(JSON.stringify(scene.models[i].animation));
                    }
                }
            } else if (model.type === 'cube') {
                this.generateModel(scene.models[i], model);
            } else if (model.type === 'cylinder') {
                this.generateModel(scene.models[i], model);
            } else if (model.type === 'cone') {
                this.generateModel(scene.models[i], model);
            } else if (model.type === 'sphere') {
                this.generateModel(scene.models[i], model);
                console.log('sphere');
            } else {
                model.center = CG.Vector4(scene.models[i].center[0],
                                       scene.models[i].center[1],
                                       scene.models[i].center[2],
                                       1);
                for (let key in scene.models[i]) {
                    if (scene.models[i].hasOwnProperty(key) && key !== 'type' && key != 'center') {
                        model[key] = JSON.parse(JSON.stringify(scene.models[i][key]));
                    }
                }
            }

            model.matrix = new Matrix(4, 4);
            processed.models.push(model);
        }
        return processed;
    }
    
    // x0:           float (x coordinate of p0)
    // y0:           float (y coordinate of p0)
    // x1:           float (x coordinate of p1)
    // y1:           float (y coordinate of p1)
    drawLine(x0, y0, x1, y1) {
        this.ctx.strokeStyle = '#000000';
        this.ctx.beginPath();
        this.ctx.moveTo(x0, y0);
        this.ctx.lineTo(x1, y1);
        this.ctx.stroke();

        this.ctx.fillStyle = '#FF0000';
        this.ctx.fillRect(x0 - 2, y0 - 2, 4, 4);
        this.ctx.fillRect(x1 - 2, y1 - 2, 4, 4);
    }
};

export { Renderer };
