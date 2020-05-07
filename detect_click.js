//Function from https://stackoverflow.com/questions/32334147/creating-hyperlinked-text-that-returns-the-text-clicked-as-input-inside-the-app
function detect_click(el) {
  Shiny.onInputChange('clicked_text', el.innerHTML);
}