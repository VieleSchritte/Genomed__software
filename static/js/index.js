/*const STEP_CHOOSE_COGNATION = 'choose_cognation';
const STEP_COGNATION_PREVIEW = 'cognation_preview';
const STEP_FILL_IN_DATA = 'fill_in_data';

let value = 0;
let currentFormula = '';

let step = STEP_CHOOSE_COGNATION;
let fillInPart = 0;
let data = [];
let participantsNumber = 2;

function HideItem(args) {
    let idList = args[0];
    for (let id of idList) {
        let item = document.getElementById(id);
        item.style.display = 'none';

    }
    if (args.length > 1) {
        let numbers = args[1];
        for (let number of numbers) {
            let id = 'part' + String(number);
            let divContainsTextarea = document.getElementById(id);
            let genotypeHeader = divContainsTextarea.querySelector('.genotype-header');
            genotypeHeader.innerHTML = ' ';
            divContainsTextarea.style.display = 'none';
        }
    }
}

function ShowItem(args) {
    let idList = args[0];
    for (let id of idList) {
        let item = document.getElementById(id);
        item.style.display = 'block';
        if (id === 'next-button-genotype' || id === 'submit-button') {
            item.style.display = 'inline-block';
        }
    }

    if (args.length > 1) {
        let numbers = args[1];
        for (let number of numbers) {
            let id = 'part' + String(number);
            let divContainsTextarea = document.getElementById(id);
            let genotypeHeader = divContainsTextarea.querySelector('.genotype-header');

            divContainsTextarea.style.display = 'inline-block';
            genotypeHeader.innerHTML = args[2][number - 1];
        }
    }
}

function GetFirstSecondGenotypes(participantsNumber) {
    if (participantsNumber % 2 === 0) {
        return {'first': [1, 2], 'second': [3, 4]};
    } else {
        return {'first': [1, 2, 3], 'second': [4, 5]};
    }
}

function ClearTextAreaValue() {
    for (let textarea of document.querySelectorAll('textarea.form-control')) {
        textarea.value = "";
    }
}

function OnRadioClick(event) {
    if (this !== event.target) {
        return
    }

    step = STEP_COGNATION_PREVIEW;
    let label = event.target;
    let input = label.querySelector('input');
    value = Number(input.getAttribute('value'));
    let imageUrl = input.getAttribute('data-image-url');
    currentFormula = formulas.find((f) => f['value'] === value);

    if (currentFormula['headers'].length > 2) {
        participantsNumber = currentFormula['headers'].length;
    }

    let cognationPreview = document.getElementById('cognation-preview');
    cognationPreview.querySelector('.preview-header').innerHTML = currentFormula['name'];
    let imageForm = cognationPreview.querySelector('.preview-image');
    imageForm.setAttribute('src', imageUrl);
    imageForm.setAttribute('alt', currentFormula['name']);
    let previewText = cognationPreview.querySelector('.preview-text');
    previewText.innerHTML = currentFormula['desc'];
    HideItem([['cognation-screen', 'start-header']]);
    ShowItem([['cognation-preview']]);
}


window.onload = function setUpHandlers() {
    for (let label of document.querySelectorAll('label.cognation-button')) {
        label.addEventListener("click", OnRadioClick);
        ClearTextAreaValue();
    }
}

function OnBackClick(event) {
    if (this !== event.target) {
        return;
    }

    if (step === STEP_COGNATION_PREVIEW) {
        step = STEP_CHOOSE_COGNATION;
        currentFormula = '';
        HideItem([['cognation-preview']]);
        ShowItem([['cognation-screen', 'start-header']]);
        ClearTextAreaValue();
        value = 0;
        currentFormula = '';
        fillInPart = 0;
        participantsNumber = 2;
        return;
    }

    if (step === STEP_FILL_IN_DATA) {
        if (fillInPart === 1) {
            step = STEP_COGNATION_PREVIEW;
            fillInPart--;
            HideItem([['add-genotype'], GetFirstSecondGenotypes(participantsNumber)['first'],
                currentFormula['headers']]);
            ShowItem([['cognation-preview']]);
        }

        if (fillInPart === 2) {
            fillInPart--;
            HideItem([['submit-button'], GetFirstSecondGenotypes(participantsNumber)['second']]);
            ShowItem([['next-button-genotype'], GetFirstSecondGenotypes(participantsNumber)['first'],
                currentFormula['headers']])
        }
    }
}

let previewBackButton = document.getElementById('back-button-preview');
previewBackButton.addEventListener('click', OnBackClick);
let backButtonGenotype = document.getElementById('back-button-genotype');
backButtonGenotype.addEventListener('click', OnBackClick);

function OnNextClick(event) {
    if (this !== event.target) {
        return
    }

    fillInPart++;

    if (step === STEP_COGNATION_PREVIEW) {
        step = STEP_FILL_IN_DATA;
        if (participantsNumber < 4) {
            HideItem([['next-button-genotype', 'cognation-preview']]);
            ShowItem([['add-genotype', 'submit-button'], GetFirstSecondGenotypes(participantsNumber)['first'],
                currentFormula['headers']]);
            return;
        } else {
            HideItem([['submit-button', 'cognation-preview']]);
            ShowItem([['add-genotype', 'next-button-genotype'], GetFirstSecondGenotypes(participantsNumber)['first'],
                currentFormula['headers']]);
            return;
        }
    }

    if (step === STEP_FILL_IN_DATA) {
        HideItem([['next-button-genotype'], GetFirstSecondGenotypes(participantsNumber)['first']]);
        ShowItem([['submit-button'], GetFirstSecondGenotypes(participantsNumber)['second'],
            currentFormula['headers']]);
    }
}

let previewNextButton = document.getElementById('next-button-preview');
previewNextButton.addEventListener("click", OnNextClick);
let nextButtonGenotype = document.getElementById('next-button-genotype');
nextButtonGenotype.addEventListener("click", OnNextClick);
*/