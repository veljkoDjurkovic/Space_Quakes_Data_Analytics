function toggleContent(contentId) {
    const contentDiv = document.getElementById(contentId);
    const allContentDivs = document.querySelectorAll('.content');

    const isAlreadyVisible = !contentDiv.classList.contains('hidden');

    if (!isAlreadyVisible) {
        allContentDivs.forEach(div => {
            div.classList.add('hidden');
        });
    }

    if (isAlreadyVisible) {
        contentDiv.classList.add('hidden'); 
    } else {
        contentDiv.classList.remove('hidden'); 
    }
}
