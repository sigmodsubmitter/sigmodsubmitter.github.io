// The following is for the show-more show-less buttons
$(".show-more a").each(function() {
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



$(".show-more a").on("click", function() {
    var $link = $(this);
    var $content = $link.parent().prev("span.text-content");
    var linkText = $link.text();
    var id = $(this)[0].id; 
    $content.toggleClass("short-text, full-text");

    $link.text(getShowLinkText(linkText, id));

    return false;
});



function getShowLinkText(currentText, id) {
    var newText = '';
    var caption;
    if (id == "caption1") {
       caption=document.getElementById("figure1_caption");
       if (currentText.toUpperCase() === "SHOW MORE...") {
          newText = "Show less.";
          caption.innerHTML = "<strong>Figure 1.</strong> The LSM-tree design space exhibits a trade-off between lookup cost and update cost that can be navigated by tuning the merge policy and size ratio. In general, the curve for Monkey  dominates the curve for the state-of-the-art because Monkey minimizes worst-case query cost by allocating main memory among the Bloom filters so as to minimize the sum of their false positive rates. To generate these curves, we plug different combinations of the merge policy and the size ratio into our closed-form equations for lookup cost and update cost, and we plot lookup cost against update cost for corresponding values of the merge policy and size ratio.  ";
      } else {
          newText = "Show more...";
          caption.innerHTML = "<strong>Figure 1.</strong> The LSM-tree design space exhibits a trade-off between lookup cost and update cost that can be navigated by tuning the merge policy and size ratio. ";
      }
    }
    else  {
       caption=document.getElementById("figure2_caption");
       
          if (currentText.toUpperCase() === "SHOW MORE...") {
          newText = "Show less.";
          caption.innerHTML = "<strong>Figure 2.</strong> Monkey minimizes lookup cost by allocating memory to Bloom filters across different levels so as to minimize the sum of their false positive rates. On the right-hand side of the figure, we display the number of entries in each level in the LSM-tree resulting from the chosen configuration. On the left-hand side, we compare how state-of-the-art designs and Monkey set the false positive rates for Bloom filters across different levels. ";
       } else {
          newText = "Show more...";
          caption.innerHTML = "<strong>Figure 2.</strong> Monkey minimizes lookup cost by allocating memory to Bloom filters across different levels so as to minimize the sum of their false positive rates.  ";
      }
    }

    return newText;
}