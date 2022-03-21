import Module from "./gauss.cjs";

const gauss = await Module();

function scopeCreate() {
    return {
        context: [],
    };
}

function scopeDestroy(scope) {
	  for(let a of scope.context) {
        a.delete();
    }
}

function number(scope, v) {
		let t = gauss.numberFromDouble(v);
		scope.context.push(t);
		return t;
}

function symbol(scope, x) {
		let t = gauss.symbol(x);
		scope.context.push(t);
		return t;
}

function add(scope, a, b) {
    let t = gauss.add(a, b);
    scope.context.push(t);
    return t;
}

function sub(scope, a, b) {
    let t = gauss.sub(a, b);
    scope.context.push(t);
    return t;
}

function mul(scope, a, b) {
    let t = gauss.mul(a, b);
    scope.context.push(t);
    return t;
}

function mul(scope, a, b) {
    let t = gauss.mul(a, b);
    scope.context.push(t);
    return t;
}

function root(scope, a, b) {
    let t = gauss.root(a, b);
    scope.context.push(t);
    return t;
}

function sqrt(scope, a) {
    let t = gauss.sqrt(a);
    scope.context.push(t);
    return t;
}

function getOperand(scope, a, i) {
		let t = gauss.getOperand(a, i);
		scope.context.push(t);
		return t;
}

function setOperand(a, i, b) {
		return gauss.setOperand(a, i, b);
}

function kindOf(a) {
		return gauss.kindOf(a);
}

function is(a, k) {
		return gauss.is(a, k);
}

function reduce(scope, a) {
		let t = gauss.reduce(a);
		scope.context.push(t);
		return t;
}

function expand(scope, a) {
		let t = gauss.expand(a);
		scope.context.push(t);
		return t;
}

function log(scope, a, b) {
    let t = gauss.log(a, b);
    scope.context.push(t);
    return t;
}

function ln(scope, a) {
    let t = gauss.ln(a);
    scope.context.push(t);
    return t;
}

function exp(scope, a) {
    let t = gauss.exp(a);
    scope.context.push(t);
    return t;
}

function replace(scope, a, x, v) {
		let t = gauss.replace(a, x, v);
		scope.context.push(t);
		return t;
}

function eval(scope, a, x, v) {
		let t = gauss.eval(a, x, v);
		scope.context.push(t);
		return t;
}

function freeVariables(scope, u) {
		let t = gauss.freeVariables(u);
		scope.context.push(t);
		return t;
}

function sinh(scope, x) {
    let t = gauss.sinh(x);
    scope.context.push(t);
    return t;
}

function cosh(scope, x) {
    let t = gauss.sinh(x);
    scope.context.push(t);
    return t;
}

function tanh(scope, x) {
    let t = gauss.tanh(x);
    scope.context.push(t);
    return t;
}

function cos(scope, x) {
    let t = gauss.cos(x);
    scope.context.push(t);
    return t;
}

function sin(scope, x) {
    let t = gauss.sin(x);
    scope.context.push(t);
    return t;
}

function tan(scope, x) {
    let t = gauss.tan(x);
    scope.context.push(t);
    return t;
}

function csc(scope, x) {
    let t = gauss.csc(x);
    scope.context.push(t);
    return t;
}

function cot(scope, x) {
    let t = gauss.cot(x);
    scope.context.push(t);
    return t;
}

function sec(scope, x) {
    let t = gauss.sec(x);
    scope.context.push(t);
    return t;
}

function coth(scope, x) {
    let t = gauss.coth(x);
    scope.context.push(t);
    return t;
}

function sech(scope, x) {
    let t = gauss.sech(x);
    scope.context.push(t);
    return t;
}

function csch(scope, x) {
    let t = gauss.csch(x);
    scope.context.push(t);
    return t;
}

function arccos(scope, x) {
    let t = gauss.arccos(x);
    scope.context.push(t);
    return t;
}

function arcsin(scope, x) {
    let t = gauss.arcsin(x);
    scope.context.push(t);
    return t;
}

function arctan(scope, x) {
    let t = gauss.arctan(x);
    scope.context.push(t);
    return t;
}

function arccot(scope, x) {
    let t = gauss.arccot(x);
    scope.context.push(t);
    return t;
}

function arcsec(scope, x) {
    let t = gauss.arcsec(x);
    scope.context.push(t);
    return t;
}

function arccsc(scope, x) {
    let t = gauss.arccsc(x);
    scope.context.push(t);
    return t;
}

function arccosh(scope, x) {
    let t = gauss.arccosh(x);
    scope.context.push(t);
    return t;
}

function arctanh(scope, x) {
    let t = gauss.arctanh(x);
    scope.context.push(t);
    return t;
}

function matrixZeros(scope, rows, columns) {
		let t = gauss.matrix(rows, columns);
		scope.context.push(t);
		return t;
}

function matrixIndentity(scope, rows, columns) {
		let t = gauss.identity(rows, columns);
		scope.context.push(t);
		return t;
}

function matrixGet(scope, M, row, column) {
		let t = gauss.matrixGet(M, row, column);
		scope.context.push(t);
		return t;
}

function matrixSet(scope, M, row, column, a) {
		let t = gauss.matrixSet(M, row, column, a);
		scope.context.push(t);
		return t;
}

function matrixSVD(scope, M) {
		let t = gauss.svd(M);
		scope.context.push(t);
		return t;
}

function matrixInverse(scope, M) {
		let t = gauss.inverse(M);
		scope.context.push(t);
		return t;
}

function matrixDeterminant(scope, M) {
		let t = gauss.det(M);
		scope.context.push(t);
		return t;
}

function matrixTranspose(scope, M) {
		let t = gauss.transpose(M);
		scope.context.push(t);
		return t;
}

function solveLinearSystem(scope, M, b) {
		let t = gauss.solveLinear(M, b);
		scope.context.push(t);
		return t;
}

function polynomialDegree(scope, p, x) {
		let t = gauss.degreePoly(p, x);
		scope.context.push(t);
		return t;
}

function polynomialCoefficient(scope, p, x, d) {
		let t = gauss.coefficientPoly(p, x, d);
		scope.context.push(t);
		return t;
}

function polynomialLeadingCoefficient(scope, p, x) {
		let t = gauss.leadingCoefficientPoly(p, x);
		scope.context.push(t);
		return t;
}

function polynomialRoots(scope, p) {
		let t = gauss.rootsOfPoly(p, x);
		scope.context.push(t);
		return t;
}


function polynomialFactors(scope, p) {
		let t = gauss.factorPoly(p);
		scope.context.push(t);
		return t;
}

function polynomialResultant(scope, u, v) {
		let t = gauss.resultantOfPoly(u, v);
		scope.context.push(t);
		return t;
}

function polynomialDivision(scope, u, v) {
		let t = gauss.divPoly(u, v);
		scope.context.push(t);
		return t;
}

function polynomialGCD(scope, u, v) {
		let t = gauss.gcdPoly(u, v);
		scope.context.push(t);
		return t;
}

function polynomialLCM(scope, u, v) {
		let t = gauss.lcmPoly(u, v);
		scope.context.push(t);
		return t;
}

function polynomialProjectFiniteField(scope, u, p) {
		let t = gauss.projectPolyFiniteField(u, u);
		scope.context.push(t);
		return t;
}

function polynomiaDivisionFiniteField(scope, u, v, p) {
		let t = gauss.projectPolyFiniteField(u, u, v);
		scope.context.push(t);
		return t;
}

function diff(scope, u, x) {
		let t = gauss.derivative(u, x);
		scope.context.push(t);
		return t;
}

function toString(a) {
		return gauss.toString(a);
}

function toLatex(a) {
		return gauss.toLatex(a);
}

function prime(scope, i) {
		let t = gauss.prime(i);
		scope.context.push(t);
		return t;
}

function primeFactors(scope, a) {
		let t = gauss.primeFactors(a);
		scope.context.push(t);
		return t;
}

let scope = scopeCreate();

let a = number(scope, 4);

let b = number(scope, 0.5);

let c = number(scope, 0.6);

let d = number(scope, 0.333333333333);

let e = symbol(scope, "x");

let f = add(scope, a, b);

let g = add(scope, f, e);

let h = add(scope, g, d);

console.log(gauss.toString(h));

scopeDestroy(scope);

gauss.doLeakCheck();
