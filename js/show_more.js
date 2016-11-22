// The following is for the show-more show-less buttons
$(".show-more_fig1 a").each(function() {
    var $link = $(this);
    var $content = $link.parent().prev("span.text-content");

    console.log($link);

    var visibleHeight = $content[0].clientHeight;
    var actualHide = $content[0].scrollHeight - 1;

    console.log(actualHide);
    console.log(visibleHeight);

    if (actualHide > visibleHeight) {
        $link.show();
    } else {
        // $link.hide();
    }
});

$(".show-more_fig1 a").on("click", function() {
    var $link = $(this);
    var $content = $link.parent().prev("span.text-content");
    var linkText = $link.text();

    $content.toggleClass("short-text, full-text");

    $link.text(getShowLinkText(linkText));

    return false;
});

function getShowLinkText(currentText) {
    var newText = '';
    var caption=document.getElementById("figure1_caption");
    if (currentText.toUpperCase() === "SHOW MORE...") {
        newText = "Show less.";
        caption.innerHTML = "<strong>Figure 1.</strong> The LSM-tree design space exhibits a trade-off between lookup cost and update cost that can be navigated by tuning the merge policy and size ratio. In general, the curve for Monkey  dominates the curve for the state-of-the-art because Monkey minimizes worst-case query cost by allocating main memory among the Bloom filters so as to minimize the sum of their false positive rates. To generate these curves, we plug different combinations of the merge policy and the size ratio into our closed-form equations for lookup cost and update cost, and we plot lookup cost against update cost for corresponding values of the merge policy and size ratio.  ";
    } else {
        newText = "Show more...";
        caption.innerHTML = "<strong>Figure 1.</strong> The LSM-tree design space exhibits a trade-off between lookup cost and update cost that can be navigated by tuning the merge policy and size ratio. ";


    }

    return newText;
}