'use strict'

function show() {
    let allDivs = document.getElementsByClassName('content');
    for (let i = 0; i < allDivs.length; i++) {
        allDivs[i].hidden = allDivs[i].hidden != true;
    }
    //document.getElementById('content' + id).classList.remove('hidden');
}