const STEP_CHOOSE_COGNATION = 'choose_cognation';
const STEP_COGNATION_PREVIEW = 'cognation_preview';
const STEP_FILL_IN_DATA = 'fill_in_data';

let step = STEP_CHOOSE_COGNATION;
/* variable for comparing data.length and participants number (is data list full or not) */
let fillInPart = 0;
let data = [];

function showCognationScreen() {
    let startHeader = document.getElementById('start-header')
    let cognationScreen = document.getElementById('cognation-screen');
    cognationScreen.style.display = 'block';
    startHeader.style.display = 'block';
}

function hideChooseCognationScreen() {
    let startHeader = document.getElementById('start-header');
    let cognationScreen = document.getElementById('cognation-screen');
    cognationScreen.style.display = 'none';
    startHeader.style.display = 'none';
}

function hideCognationPreviewScreen() {
    let cognationPreviewScreen = document.getElementById('cognation-preview')
    cognationPreviewScreen.style.display = 'none';
}

function showAddGenotypeScreen() {
    let addGenotypeScreen = document.getElementById('add-genotype');
    addGenotypeScreen.style.display = 'block';
}

function hideAddGenotypeScreen() {
    let addGenotypeScreen = document.getElementById('add-genotype');
    addGenotypeScreen.style.display = 'none';
}

function showSubmitButton() {
    let submitButton = document.getElementById('submit-button');
    submitButton.style.display = 'inline-block';
}

function hideSubmitButton() {
    let submitButton = document.getElementById('submit-button');
    submitButton.style.display = 'none';
}


function onBackClick(event) {
    if (step === STEP_COGNATION_PREVIEW) {
        step = STEP_CHOOSE_COGNATION;
        /* hide cognation preview, show radio buttons */
        hideCognationPreviewScreen();
        showCognationScreen();
    }
    if (step === STEP_FILL_IN_DATA) {
        if (fillInPart === 0) {
            /* on the first fill in => show preview, hide fill in */
            step = STEP_CHOOSE_COGNATION;
            showCognationPreviewScreen();
            hideAddGenotypeScreen();
        }
        if (fillInPart < formulas['participants_number']) {
            /* on other fill ins => show previous fill in, hide current */
            fillInPart --;
        }
    }
}

function onNextClick() {
    if (step === STEP_COGNATION_PREVIEW) {
        step = STEP_FILL_IN_DATA;
        /* hide cognation preview, show fill in first */
        hideCognationPreviewScreen;
        showAddGenotypeScreen;
    }
    if (step === STEP_FILL_IN_DATA) {
        fillInPart ++;
        /* hide current fill in, show next */
    }
}

function onSubmitClick(event) {

}

/* Behaviour of radio buttons - on click on one of them show its description and hide the rest of buttons */
for (let label of document.querySelectorAll('label.cognation-button')) {
    step = STEP_COGNATION_PREVIEW;
    let labelsValue = label.firstElementChild.getAttribute('value')
    label.onclick = hideChooseCognationScreen(labelsValue)
}

/* On click on back button - hide description and show the rest of radio buttons */
for (let backButton of document.querySelectorAll('label.description-back')) {
    backButton.onclick = function OnBackButtonClick() {
        backButton.closest('.description-genotypes').style.display = 'none';

        for (let label of document.querySelectorAll('label.cognation-button')) {
            label.style.display = 'inline';
            label.classList.remove("active");
            let startHeader = document.getElementById('start-header');
            startHeader.style.display = 'block';
        }
    }
}

/* On click on ok button in a description - hide description and label, show first field fof genotype, hide start header */
for (let okButton of document.querySelectorAll('label.description-ok')) {
    okButton.onclick = function OnOKButtonClick() {
        /* hiding label */
        let labelToHide = okButton.closest('.description-genotypes').previousElementSibling;
        labelToHide.style.display = 'none';

        /* hiding description */
        let cognDescr = okButton.closest('.description')
        cognDescr.style.display = 'none';

        /* hiding start header */
        document.getElementById('start-header').style.display = 'none';

        /* showing first field for genotype */
        let addGenotype = cognDescr.nextElementSibling.firstElementChild;
        addGenotype.style.display = 'block';
    }
}

/* First child back - hides itself, shows description */
for (let allGenotypes of document.querySelectorAll('div.all-genotypes')) {
    let buttonGroup = allGenotypes.querySelector('div.first-genotype').querySelector('div.genotype-button-group');
    let firstBackButton = buttonGroup.querySelector('label.first-back');
    firstBackButton.onclick = function firstBackButtonBehaviour() {
        let firstGenotype = allGenotypes.firstElementChild;
        firstGenotype.style.display = 'none';

        let description = allGenotypes.previousElementSibling;
        description.style.display = 'block';
    }
}

/* All other back buttons - hide themselves, show previous genotypes */
for (let genotypeBackButton of document.querySelectorAll('label.genotype-back')) {
    let buttonGroup = genotypeBackButton.closest('.genotype-button-group');
    let currentGenotype = buttonGroup.parentElement;
    let previousGenotype = currentGenotype.previousElementSibling;
    genotypeBackButton.onclick = function OnBackClick() {
        currentGenotype.style.display = 'none';
        previousGenotype.style.display = 'block';
    }
}

/* Next button always shows next genotype */
for (let genotypeNextButton of document.querySelectorAll('label.genotype-next')) {
                let buttonGroup = genotypeNextButton.closest('.genotype-button-group');
                let currentGenotype = buttonGroup.parentElement;
                let nextGenotype = currentGenotype.nextElementSibling;
                genotypeNextButton.onclick = function OnNextClick() {
                    currentGenotype.style.display = 'none';
                    nextGenotype.style.display = 'block';
                }
            }