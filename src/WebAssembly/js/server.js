// import * as http from 'http'; //ES 6
// import * as fs from 'fs'; //ES 6

const host = 'localhost';
const port = 8000;

const http = require("http");

const Module = require("./gaussjs.js");

Module().then(function(instance) {
		console.log(instance)
});

//let t = require('./gaussjs.js');
//console.log(t.exports)
// import * as k from './gaussjs.wasm.js'
        // var Module = {
        //     onRuntimeInitialized: function() {
        //         let n = Module.num(5)

        //         console.log('lerp result: ' + Module.toString(n));

        //         n.delete()

        //     }
        // };



const requestListener = function (req, res) {
//     console.log("request:" + req.url)

//     switch (req.url) {
//         case "/":
//             fs.readFile(__dirname + "/index.html")
//                 .then(contents => {
//                     res.setHeader("Content-Type", "text/html");
//                     res.writeHead(200);
//                     res.end(contents);
//                 })
//                 .catch(err => {
//                     res.writeHead(500);
//                     res.end(err);
//                     return;
//                 });
//             break
//         case "/gaussjs.wasm.js":
//             fs.readFile(__dirname + "/gaussjs.wasm.js")
//             .then(contents => {
//                 res.setHeader("Content-Type", "text/javascript");
//                 res.writeHead(200);
//                 res.end(contents);
//             })
//             .catch(err => {
//                 res.writeHead(500);
//                 res.end(err);
//                 return;
//             });
//         break
//         case "/gaussjs.wasm.wasm":
//             fs.readFile(__dirname + "/gaussjs.wasm.wasm")
//             .then(contents => {
//                 res.setHeader("Content-Type", "application/wasm");
//                 res.writeHead(200);
//                 res.end(contents);
//             })
//             .catch(err => {
//                 res.writeHead(500);
//                 res.end(err);
//                 return;
//             });
//         break
//         case "/gauss.js":
//             fs.readFile(__dirname + "/gauss.js")
//             .then(contents => {
//                 res.setHeader("Content-Type", "text/javascript");
//                 res.writeHead(200);
//                 res.end(contents);
//             })
//             .catch(err => {
//                 res.writeHead(500);
//                 res.end(err);
//                 return;
//             });
//         break
//         default:
//             res.writeHead(404);
//             res.end(JSON.stringify({error:"Resource not found"}));
//     }
}

const server = http.createServer(requestListener);

server.listen(port, host, () => {
     console.log(`Server is running on http://${host}:${port}`);
});
