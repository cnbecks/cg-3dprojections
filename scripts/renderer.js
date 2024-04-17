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
        this.enable_animation = true;  // <-- disabled for easier debugging; enable for animation
        this.start_time = null;
        this.prev_time = null;
        this.rotated = false;
    }

    //
    updateTransforms(time, delta_time) {
        // using time rather than delta time...
        // NOT dealing with a list of all models, one model at a time
        let model = this.scene.models[0];
        if ( model.hasOwnProperty('animation') ) {
            // set transformation matrix to transport model to the origin
            let transform = model.animation.transform;
            transform = new Matrix(4, 4);
            CG.mat4x4Translate( transform, -1*model.center[0], -1*model.center[1], -1*model.center[2] );

            // calculate rotation matrix based on model velocity in rps
            let theta = ( model.animation.rps / 1000.0 ) * delta_time; // time is in milliseconds
            if ( theta >= (2*Math.PI) ) {
                theta = theta - (2*Math.PI); // resets theta after a full circle is complete
            }
            // determine which axis the rotation is in to complete matrix calculation
            let rot_matrix = new Matrix(4, 4);
            if ( model.animation.axis == "x" ) {
                CG.mat4x4RotateX( rot_matrix, theta );
            } else if ( model.animation.axis == "y" ) {
                CG.mat4x4RotateY( rot_matrix, theta );
            } else { // model.animation.axis == "z"
                CG.mat4x4RotateZ( rot_matrix, theta );
            }

            //loop through each vertex in the model and apply the series of transforms to each one
            let updated_vertices = [];
            for (let j = 0; j < model.vertices.length; j++) { // 0-4 = back, 5-9 = front
                // calculate reverse translation matrix
                let rev_transform = new Matrix(4, 4);
                CG.mat4x4Translate( rev_transform, model.center[0], model.center[1], model.center[2] );

                updated_vertices.push( Matrix.multiply([rev_transform, rot_matrix, transform, model.vertices[j]]) )    
            }
            model.vertices = updated_vertices;
        }  
    }

    generateModel(model, new_model) {
        let center = model.center;
        let vertices = [];
        let edges = [];

        //---------------- CUBE ----------------//
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

        //---------------- CONE ----------------//
        } else if (model.type == "cone") {
            let h = model.height/2;
            let r = model.radius;
            
            // generate circular face
            for (let i = 0; i < model.sides; i++) {
                let currentTheta = (i / model.sides) * (2 * Math.PI);
                vertices.push( CG.Vector4( model.center[0] + model.radius * Math.cos(currentTheta),
                                           model.center[1] + h,
                                           model.center[2] + model.radius * Math.sin(currentTheta), 
                                           1 ) );
            }

            // push cone tip coordinate
            vertices.push( CG.Vector4( model.center[0],
                                       model.center[1] - h,
                                       model.center[2], 
                                       1 ) );

            // create edge matrices
            let cir_edge = [];
            for (let num_sd = 0; num_sd < model.sides; num_sd++) {
                cir_edge.push(num_sd);                         // create circle face edge
                edges.push( [cir_edge[num_sd], model.sides] ); // create top to bottom edges
            }
            cir_edge.push(0);
            edges.push(cir_edge);

        //---------------- CYLINDER ----------------//
        } else if (model.type == "cylinder") {
            let h = model.height/2;
            let r = model.radius;
            
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

            // create edge matrices
            let top_edge = [];
            let bot_edge = [];
            for (let num_sd = 0; num_sd < model.sides; num_sd++) {
                top_edge.push(num_sd);                              // create top circle edges
                bot_edge.push(model.sides + num_sd);                // create bottom circle edges
                edges.push( [top_edge[num_sd], bot_edge[num_sd]] ); // create top to bottom edges
            }
            top_edge.push(0);
            edges.push(top_edge);
            bot_edge.push(model.sides);
            edges.push(bot_edge);

        //---------------- SPHERE ----------------//
        } else if (model.type == "sphere") {
            let r = model.radius;
            let h = (2 * r) / model.stacks;
            
            // generate top and bottom circular faces
            for (let stack = 1; stack < model.stacks; stack++) {
                let stack_edge = [];
                let middle = Math.floor( model.stacks / 2 );
                for (let i = 0; i < model.slices; i++) {
                    let current_radius;
                    if ( stack == middle && (model.stacks / 2) != middle) { // middle of the sphere with an ODD number of stacks
                        current_radius = r;
                    } else if ( stack > middle || (stack == middle && (model.stacks / 2) == middle) ) { // top half of the sphere
                        current_radius = Math.round( Math.sqrt( (r**2) - ((h*stack)-r)**2 ) ); // solve for b in c^2 = a^2 + b^2
                    } else { // bottom half of the sphere
                        current_radius = Math.round( Math.sqrt( (r**2) - (r-(h*stack))**2 ) ); // solve for b in c^2 = a^2 + b^2
                    }

                    let currentTheta = (i / (model.slices)) * (2 * Math.PI);
                    vertices.push( CG.Vector4( model.center[0] + current_radius * Math.cos(currentTheta),
                                               model.center[1] - r + (h*stack),
                                               model.center[2] + current_radius * Math.sin(currentTheta), 
                                               1 ) );
                    
                    // add edges for this stack
                    stack_edge.push( (model.slices*stack) + i - model.slices ); // add 2 to take into account top and bottom vertices
                }
                stack_edge.push(stack_edge[0]);
                // add stack edges
                edges.push(stack_edge);
            }

            // push top and bottom vertices
            vertices.push( CG.Vector4( model.center[0],
                                       model.center[1] - r,
                                       model.center[2], 
                                       1 ) );

            vertices.push( CG.Vector4( model.center[0],
                                       model.center[1] + r,
                                       model.center[2], 
                                       1 ) );

            // get top and nottom edge positions
            let top_loc = vertices.length-1;
            let bot_loc = vertices.length-2;

            // add slice edges
            for (let slice = 0; slice < model.slices-1; slice++) {
                let slice_edge = [];
                slice_edge.push(bot_loc);
                for (let i = 0; i < model.stacks-1; i++) {
                    slice_edge.push(edges[i][slice]);
                }
                slice_edge.push(top_loc);
                edges.push(slice_edge);
            }
        }
        new_model.vertices = vertices;
        new_model.edges = edges;
        new_model.center = center;
    }
    

    rotateLeft() {
        let omega = 10 * Math.PI/180;
        let prp = this.scene.view.prp;
        let srp = this.scene.view.srp;

        // translate SRP to origin/prp
        let translation_matrix = new Matrix(4,4);
        CG.mat4x4Translate(translation_matrix, -prp.x, -prp.y, -prp.z); // set translation matrix to move the srp back to the origin

        // rotate VRC such that (u,v,n) align with (x,y,z)
        let n = prp.subtract(srp);
        n.normalize();
        let u = this.scene.view.vup.cross(n);
        u.normalize();
        let v = n.cross(u);
        let R = new Matrix(4,4);
        R.values = [[u.x, u.y, u.z, 0],
                [v.x, v.y, v.z, 0],
                [n.x, n.y, n.z, 0],
                [0,    0,    0,    1]];

        let align = Matrix.multiply([R, translation_matrix]);

        let rotation_y = new Matrix(4,4);
        CG.mat4x4RotateY(rotation_y, omega); //rotate it based on the y-axis now that it is aligned with the v axis

        // ------------ UNDO THE VRC rotation and the translation_matrix ------------------------------
        // undo the the VRC alignement - By setting these values to the opposite of what they were before
        // let undo_R = R.inverse();
        // let undo_R = R.inverse();

        // undo the translation_matrix
        // let undo_translation_matrix = translation_matrix.inverse(); //new Matrix(4,4);
        let undo = align.inverse(); //new Matrix(4,4);
        
        // now apply these transformations to the srp
        let vector_4 = CG.Vector4(srp.x, srp.y, srp.z, 1); //vector_4 will hold the new values of the srp
        let translations = Matrix.multiply([undo, rotation_y, align, vector_4]);

        // Set the SRP to the new translated, rotated, then translated back to the srp 
        this.scene.view.srp.x = translations.x/translations.w;
        this.scene.view.srp.y = translations.y/translations.w;
        this.scene.view.srp.z = translations.z/translations.w;

        // hack Dr. Marrinan showed me - call rotateLeft TWICE to avoid the house disappearing problem
        if (this.rotated == false) { //this.rotated is a global variable
            this.rotated = true;
            this.rotateLeft();
        }
        this.rotated = false;
        this.draw();
    }
    
    

    rotateRight() {
        let omega = -10 * Math.PI/180;
        let prp = this.scene.view.prp;
        let srp = this.scene.view.srp;

        // translate SRP to origin/prp
        let translation_matrix = new Matrix(4,4);
        CG.mat4x4Translate(translation_matrix, -prp.x, -prp.y, -prp.z); // set translation matrix to move the srp back to the origin

        // rotate VRC such that (u,v,n) align with (x,y,z)
        let n = prp.subtract(srp);
        n.normalize();
        let u = this.scene.view.vup.cross(n);
        u.normalize();
        let v = n.cross(u);
        let R = new Matrix(4,4);
        R.values = [[u.x, u.y, u.z, 0],
                [v.x, v.y, v.z, 0],
                [n.x, n.y, n.z, 0],
                [0,    0,    0,    1]];

        let align = Matrix.multiply([R, translation_matrix]);

        let rotation_y = new Matrix(4,4);
        CG.mat4x4RotateY(rotation_y, omega); //rotate it based on the y-axis now that it is aligned with the v axis

        // ------------ UNDO THE VRC rotation and the translation_matrix ------------------------------
        // undo the the VRC alignement - By setting these values to the opposite of what they were before
        // let undo_R = R.inverse();

        // undo the translation_matrix
        // let undo_translation_matrix = translation_matrix.inverse(); //new Matrix(4,4);
        let undo = align.inverse(); //new Matrix(4,4);
        
        // now apply these transformations to the srp
        let vector_4 = CG.Vector4(srp.x, srp.y, srp.z, 1); //vector_4 will hold the new values of the srp
        let translations = Matrix.multiply([undo, rotation_y, align, vector_4]);

        // Set the SRP to the new translated, rotated, then translated back to the srp 
        this.scene.view.srp.x = translations.x/translations.w;
        this.scene.view.srp.y = translations.y/translations.w;
        this.scene.view.srp.z = translations.z/translations.w;

        // hack Dr. Marrinan showed me - call rotateLeft TWICE to avoid the house disappearing problem
        if (this.rotated == false) { //this.rotated is a global variable
            this.rotated = true;
            this.rotateRight();
        }
        this.rotated = false;
        this.draw();
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
            for (let j=0; j<this.scene.models[i].vertices.length; j++) { // 0-4 = back, 5-9 = front
                let new_vertex = Matrix.multiply([perspective_matrix, this.scene.models[i].vertices[j]]);
                vertices.push(new_vertex);       
            }

            //   * For (each line segment in) each edge
            for (let k=0; k<this.scene.models[i].edges.length; k++) {

                // loop through the vertices in each edge which have number that correspond to the vertices...
                // For each line segment (in each edge)
                let actual_vertices = []; //actual vertices for the edge
                let actual_edges = [];
                for (let v=0; v<this.scene.models[i].edges[k].length; v++) {
                    // use the vertex numbers in the edge to create a list with the actual vertices
                    let vertex_number = this.scene.models[i].edges[k][v]; //in model i, in the edge k, get the vertex number v
                    actual_edges.push(this.scene.models[i].edges[k][v]);
                    let vertex_from_number = vertices[vertex_number]; // using this number, get its corresponding Vector
                    actual_vertices.push(vertex_from_number); // push the actual vertex, a Vector, onto the list
                }

                // the actual_vertices correspond to the edge --> line
                let clipped_vertices = [];
                //loop through the actual_vertices and clip each line
                for (let i=0; i<actual_vertices.length-1; i++) {
                    let line = {pt0: actual_vertices[i], pt1: actual_vertices[i+1]}
                    let z_min = (-1 * this.scene.view.clip[4]) / this.scene.view.clip[5]; // -near/far
                    let new_line = this.clipLinePerspective(line, z_min);
                    if (new_line != null) {
                        clipped_vertices.push(new_line.pt0);
                        clipped_vertices.push(new_line.pt1);
                    }
                }   
                actual_vertices = clipped_vertices;

                // Project to 2D by multiplying by MPer
                let vertices_of_edge = []
                for (let p=0; p<actual_vertices.length; p++) {
                    let mper_vertex = Matrix.multiply([mper_matrix, actual_vertices[p]]);
                    vertices_of_edge.push(mper_vertex);
                }

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
        else if (vertex.y < (vertex.z - FLOAT_EPSILON)) {
            outcode += BOTTOM;
        }
        else if (vertex.y > (-vertex.z + FLOAT_EPSILON)) {
            outcode += TOP;
        }
        else if (vertex.z < (-1.0 - FLOAT_EPSILON)) {
            outcode += FAR;
        }
        else if (vertex.z > (z_min + FLOAT_EPSILON)) {
            outcode += NEAR;
        }
        return outcode;
    }

    calcXYZ(pt0, change_x, change_y, change_z, t){
        let x = pt0.x + t*change_x;
        let y = pt0.y + t*change_y;
        let z = pt0.z + t*change_z;
        let intersect_pt = CG.Vector3(x, y, z);
        return intersect_pt;
    }

    calctLeft(pt0, pt1) { 
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = ((-1*pt0.x) + pt0.z)/(change_x-change_z);
        return this.calcXYZ(pt0, change_x, change_y, change_z, t);
    }

    calctRight(pt0, pt1) {
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = (pt0.x + pt0.z)/((-1*change_x)-change_z);
        return this.calcXYZ(pt0, change_x, change_y, change_z, t);
    }

    calctBottom(pt0, pt1) {
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = ((-1)*pt0.y + pt0.z)/(change_y-change_z);
        return this.calcXYZ(pt0, change_x, change_y, change_z, t);
    }

    calctTop(pt0, pt1) {
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = (pt0.y + pt0.z)/((-1*change_y)-change_z);
        return this.calcXYZ(pt0, change_x, change_y, change_z, t);
    }

    calctNear(pt0, pt1, z_min) {
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = (pt0.z-z_min)/(-1*change_z);
        return this.calcXYZ(pt0, change_x, change_y, change_z, t);
    }

    calctFar(pt0, pt1) {
        let change_x = (pt1.x-pt0.x);
        let change_y = (pt1.y-pt0.y);
        let change_z = (pt1.z-pt0.z);
        let t = ((-1*pt0.z)-1)/(change_z);
        return this.calcXYZ(pt0, change_x, change_y, change_z, t);
    }    

    // Clip line - should either return a new line (with two endpoints inside view volume)
    //             or null (if line is completely outside view volume)
    // line:       object {pt0: Vector4, pt1: Vector4}
    // z_min:      float (near clipping plane in canonical view volume)
    clipLinePerspective(line, z_min) {
        let result = null;
        let p0 = CG.Vector3(line.pt0.x, line.pt0.y, line.pt0.z); 
        let p1 = CG.Vector3(line.pt1.x, line.pt1.y, line.pt1.z);
        let out0 = this.outcodePerspective(p0, z_min);
        let out1 = this.outcodePerspective(p1, z_min);

        let deviation = 0.001 // to remove lines when p0 basically == p1
        while (1) {
            p0 = CG.Vector3( parseFloat(p0.x.toFixed(3)), parseFloat(p0.y.toFixed(3)), parseFloat(p0.z.toFixed(3)) );
            p1 = CG.Vector3( parseFloat(p1.x.toFixed(3)), parseFloat(p1.y.toFixed(3)), parseFloat(p1.z.toFixed(3)) );
            out0 = this.outcodePerspective(p0, z_min); // each time we loop we have to check the outcodes!
            out1 = this.outcodePerspective(p1, z_min);
            if ((out0 | out1) == 0) {
                result ={pt0: CG.Vector4(p0.x, p0.y, p0.z, 0), pt1: CG.Vector4(p1.x, p1.y, p1.z, 0)};
                return result;
            } else if ( (Math.abs(p0.x - p1.x) < deviation) && (Math.abs(p0.y - p1.y) < deviation) && (Math.abs(p0.z - p1.z) < deviation) ) {
                return result;
            } else if ((out0 & out1) != 0) { 
                return result;
            } 

            if (out0 == LEFT) {
                p0 = this.calctLeft(p0, p1);
            } else if (out0 == RIGHT) {
                p0 = this.calctRight(p0, p1);
            } else if (out0 == BOTTOM) {
                p0 = this.calctBottom(p0, p1);
            } else if (out0 == TOP) {
                p0 = this.calctTop(p0, p1);
            } else if (out0 == FAR) {
                p0 = this.calctFar(p0, p1);
            } else if ((out0 == NEAR)) { 
                p0 = this.calctNear(p0, p1, z_min);
            }

            if (out1 == LEFT) {
                p1 = this.calctLeft(p1, p0);
            } else if (out1 == RIGHT) {
                p1 = this.calctRight(p1, p0);
            } else if (out1 == BOTTOM) {
                p1 = this.calctBottom(p1, p0);
            } else if (out1 == TOP) {
                p1 = this.calctTop(p1, p0);
            } else if (out1 == FAR) {
                p1 = this.calctFar(p1, p0);
            } else if ((out1 == NEAR)) { 
                p1 = this.calctNear(p1, p0, z_min);
            }
            console.log('clipppinngg');
        }
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
                if (scene.models[i].hasOwnProperty('animation')) {
                    model.animation = JSON.parse(JSON.stringify(scene.models[i].animation));
                }
            } else if (model.type === 'cylinder') {
                this.generateModel(scene.models[i], model);
                if (scene.models[i].hasOwnProperty('animation')) {
                    model.animation = JSON.parse(JSON.stringify(scene.models[i].animation));
                }
            } else if (model.type === 'cone') {
                this.generateModel(scene.models[i], model);
                if (scene.models[i].hasOwnProperty('animation')) {
                    model.animation = JSON.parse(JSON.stringify(scene.models[i].animation));
                }
            } else if (model.type === 'sphere') {
                this.generateModel(scene.models[i], model);
                if (scene.models[i].hasOwnProperty('animation')) {
                    model.animation = JSON.parse(JSON.stringify(scene.models[i].animation));
                }
            } else {
                model.center = CG.Vector4( scene.models[i].center[0],
                                           scene.models[i].center[1],
                                           scene.models[i].center[2],
                                           1 );
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
