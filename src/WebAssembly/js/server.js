const http = require("http");
const fs   = require('fs').promises;

const host = 'localhost';
const port = 8000;

const requestListener = function (req, res) {
    console.log("request:" + req.url)

    switch (req.url) {
        case "/":
            fs.readFile(__dirname + "/index.html")
                .then(contents => {
                    res.setHeader("Content-Type", "text/html");
                    res.writeHead(200);
                    res.end(contents);
                })
                .catch(err => {
                    res.writeHead(500);
                    res.end(err);
                    return;
                });
            break
        case "/gaussjs_wasm.js":
            fs.readFile(__dirname + "/gaussjs_wasm.js")
            .then(contents => {
                res.setHeader("Content-Type", "text/javascript");
                res.writeHead(200);
                res.end(contents);
            })
            .catch(err => {
                res.writeHead(500);
                res.end(err);
                return;
            });
        break
        case "/gaussjs.wasm.wasm":
            fs.readFile(__dirname + "/gaussjs_wasm.wasm")
            .then(contents => {
                res.setHeader("Content-Type", "application/wasm");
                res.writeHead(200);
                res.end(contents);
            })
            .catch(err => {
                res.writeHead(500);
                res.end(err);
                return;
            });
        break
        default:
            res.writeHead(404);
            res.end(JSON.stringify({error:"Resource not found"}));
    }

}

const server = http.createServer(requestListener);

server.listen(port, host, () => {
    console.log(`Server is running on http://${host}:${port}`);
});