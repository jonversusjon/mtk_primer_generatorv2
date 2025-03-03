import React from "react";
import "../../styles/PrimerAnatomy.css";

function PrimerAnatomy() {
  return (
    <div className="primer-anatomy">
      <h3>BsmBI Cutting Pattern in Primers</h3>

      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 500">
        <style>
          {`
            /* Light mode colors */
            :root {
                --bg-color: white;
                --text-color: black;
                --muted-text: #666666;
                --spacer-bg: #E8E8E8;
                --bsmbi-bg: #FFE8E8;
                --seq-bg: #E8E8FF;
                --binding-bg: #F0F0F0;
                --legend-bg: #f8f8f8;
                --legend-border: #ddd;
            }
            
            /* Dark mode colors */
            :root .dark-mode {
                --bg-color: #1a1a1a;
                --text-color: #ffffff;
                --muted-text: #999999;
                --spacer-bg: #2d2d2d;
                --bsmbi-bg: #3d2828;
                --seq-bg: #28283d;
                --binding-bg: #2d2d2d;
                --legend-bg: #2d2d2d;
                --legend-border: #404040;
            }

            .bg-rect { fill: var(--bg-color); }
            .spacer-domain { fill: var(--spacer-bg); }
            .bsmbi-domain { fill: var(--bsmbi-bg); }
            .seq-domain { fill: var(--seq-bg); }
            .binding-domain { fill: var(--binding-bg); }
            .main-text { fill: var(--text-color); }
            .muted-text { fill: var(--muted-text); }
            .legend-box { fill: var(--legend-bg); stroke: var(--legend-border); }
          `}
        </style>

        {/* Background */}
        <rect width="800" height="500" className="bg-rect" />

        {/* Title */}
        <text
          x="400"
          y="40"
          textAnchor="middle"
          fontFamily="Arial"
          fontSize="20"
          fontWeight="bold"
          className="main-text"
        >
          BsmBI Cutting Pattern in Primers
        </text>

        {/* Forward Primer */}
        <g transform="translate(100, 150)">
          {/* Domain background rectangles */}
          <rect
            x="-5"
            y="-15"
            width="45"
            height="40"
            className="spacer-domain"
            rx="3"
          />
          <rect
            x="40"
            y="-15"
            width="70"
            height="40"
            className="bsmbi-domain"
            rx="3"
          />
          <rect
            x="110"
            y="-15"
            width="80"
            height="40"
            className="seq-domain"
            rx="3"
          />
          <rect
            x="190"
            y="-15"
            width="190"
            height="40"
            className="binding-domain"
            rx="3"
          />

          {/* DNA sequence */}
          <text x="0" y="0" fontFamily="monospace" fontSize="16">
            <tspan className="muted-text" x="5">
              gaa
            </tspan>
            <tspan x="50" className="main-text">
              cgtctc
            </tspan>
            <tspan x="117" className="main-text">
              G
            </tspan>
            <tspan x="132" fontWeight="bold" className="main-text">
              AAAC
            </tspan>
            <tspan x="175" className="main-text">
              C
            </tspan>
            <tspan x="200" className="muted-text">
              tccgcgccccgcaacctcc
            </tspan>
          </text>
          <text x="0" y="20" fontFamily="monospace" fontSize="16">
            <tspan className="muted-text" x="5">
              ctt
            </tspan>
            <tspan x="50" className="main-text">
              gcagag
            </tspan>
            <tspan x="117" className="main-text">
              C
            </tspan>
            <tspan x="132" fontWeight="bold" className="main-text">
              TTTG
            </tspan>
            <tspan x="175" className="main-text">
              G
            </tspan>
            <tspan x="200" className="muted-text">
              aggcgcggggcgttggagg
            </tspan>
          </text>

          {/* Forward primer cutting arrows */}
          <path
            d="M 130,-40 L 130,-10"
            stroke="#E63946"
            strokeWidth="2"
            fill="none"
            markerEnd="url(#arrowhead1)"
            className="forward-bottom"
          />
          <path
            d="M 172,50 L 172,20"
            stroke="#E63946"
            strokeWidth="2"
            fill="none"
            markerEnd="url(#arrowhead2)"
            className="forward-top"
          />

          {/* Domain labels */}
          <text
            x="15"
            y="-25"
            fontFamily="Arial"
            fontSize="12"
            textAnchor="middle"
            className="main-text"
          >
            Spacer
          </text>
          <text
            x="85"
            y="-25"
            fontFamily="Arial"
            fontSize="12"
            textAnchor="middle"
            className="main-text"
          >
            BsmBI site
          </text>
          <text
            x="175"
            y="-25"
            fontFamily="Arial"
            fontSize="12"
            textAnchor="middle"
            className="main-text"
          >
            6nt sequence
          </text>
          <text
            x="340"
            y="-25"
            fontFamily="Arial"
            fontSize="12"
            textAnchor="middle"
            className="main-text"
          >
            Binding region
          </text>

          {/* Primer label */}
          <text
            x="-80"
            y="10"
            fontFamily="Arial"
            fontSize="14"
            fontWeight="bold"
            className="main-text"
          >
            Forward
          </text>
        </g>

        {/* Reverse Primer */}
        <g transform="translate(100, 300)">
          {/* Domain background rectangles */}
          <rect
            x="-5"
            y="-15"
            width="45"
            height="40"
            className="spacer-domain"
            rx="3"
          />
          <rect
            x="40"
            y="-15"
            width="70"
            height="40"
            className="bsmbi-domain"
            rx="3"
          />
          <rect
            x="110"
            y="-15"
            width="80"
            height="40"
            className="seq-domain"
            rx="3"
          />
          <rect
            x="190"
            y="-15"
            width="190"
            height="40"
            className="binding-domain"
            rx="3"
          />

          {/* DNA sequence */}
          <text x="0" y="0" fontFamily="monospace" fontSize="16">
            <tspan className="muted-text" x="5">
              gaa
            </tspan>
            <tspan x="50" className="main-text">
              cgtctc
            </tspan>
            <tspan x="117" className="main-text">
              A
            </tspan>
            <tspan x="132" fontWeight="bold" className="main-text">
              CAGT
            </tspan>
            <tspan x="175" className="main-text">
              G
            </tspan>
            <tspan x="200" className="muted-text">
              tctcggtgaaggcctcccag
            </tspan>
          </text>
          <text x="0" y="20" fontFamily="monospace" fontSize="16">
            <tspan className="muted-text" x="5">
              ctt
            </tspan>
            <tspan x="50" className="main-text">
              gcagag
            </tspan>
            <tspan x="117" className="main-text">
              T
            </tspan>
            <tspan x="132" fontWeight="bold" className="main-text">
              GTCA
            </tspan>
            <tspan x="175" className="main-text">
              C
            </tspan>
            <tspan x="200" className="muted-text">
              agagccacttccggagggtc
            </tspan>
          </text>

          {/* Reverse primer cutting arrows */}
          <path
            d="M 130,-40 L 130,-10"
            stroke="#E63946"
            strokeWidth="2"
            fill="none"
            markerEnd="url(#arrowhead3)"
            className="reverse-bottom"
          />
          <path
            d="M 172,50 L 172,20"
            stroke="#E63946"
            strokeWidth="2"
            fill="none"
            markerEnd="url(#arrowhead4)"
            className="reverse-top"
          />

          {/* Primer label */}
          <text
            x="-80"
            y="10"
            fontFamily="Arial"
            fontSize="14"
            fontWeight="bold"
            className="main-text"
          >
            Reverse
          </text>
        </g>

        {/* Arrow marker definitions */}
        <defs>
          <marker
            id="arrowhead1"
            markerWidth="10"
            markerHeight="7"
            refX="9"
            refY="3.5"
            orient="auto"
          >
            <polygon points="0 0, 10 3.5, 0 7" fill="#E63946" />
          </marker>
          <marker
            id="arrowhead2"
            markerWidth="10"
            markerHeight="7"
            refX="9"
            refY="3.5"
            orient="auto"
          >
            <polygon points="0 0, 10 3.5, 0 7" fill="#E63946" />
          </marker>
          <marker
            id="arrowhead3"
            markerWidth="10"
            markerHeight="7"
            refX="9"
            refY="3.5"
            orient="auto"
          >
            <polygon points="0 0, 10 3.5, 0 7" fill="#E63946" />
          </marker>
          <marker
            id="arrowhead4"
            markerWidth="10"
            markerHeight="7"
            refX="9"
            refY="3.5"
            orient="auto"
          >
            <polygon points="0 0, 10 3.5, 0 7" fill="#E63946" />
          </marker>
        </defs>

        {/* Legend */}
        <g transform="translate(100, 400)">
          <rect
            x="0"
            y="0"
            width="600"
            height="80"
            className="legend-box"
            rx="5"
          />
          <text
            x="20"
            y="30"
            fontFamily="Arial"
            fontSize="14"
            fontWeight="bold"
            className="main-text"
          >
            Primer Domains:
          </text>
          <g transform="translate(20, 50)">
            <rect
              x="0"
              y="-15"
              width="30"
              height="20"
              className="spacer-domain"
              rx="3"
            />
            <text
              x="40"
              y="0"
              fontFamily="Arial"
              fontSize="12"
              className="main-text"
            >
              Spacer
            </text>

            <rect
              x="120"
              y="-15"
              width="30"
              height="20"
              className="bsmbi-domain"
              rx="3"
            />
            <text
              x="150"
              y="0"
              fontFamily="Arial"
              fontSize="12"
              className="main-text"
            >
              BsmBI recognition site
            </text>

            <rect
              x="300"
              y="-15"
              width="30"
              height="20"
              className="seq-domain"
              rx="3"
            />
            <text
              x="280"
              y="0"
              fontFamily="Arial"
              fontSize="12"
              className="main-text"
            >
              6nt sequence (4bp overhang in bold)
            </text>

            <rect
              x="480"
              y="-15"
              width="30"
              height="20"
              className="binding-domain"
              rx="3"
            />
            <text
              x="520"
              y="0"
              fontFamily="Arial"
              fontSize="12"
              className="main-text"
            >
              Binding region
            </text>
          </g>
        </g>
      </svg>

      <div className="primer-anatomy-explanation">
        <p>
          The diagram shows the structure of Golden Gate primers using BsmBI
          restriction enzyme:
        </p>
        <ul>
          <li>
            <strong>Spacer (gray):</strong> Short sequence at the 5' end to
            provide space for the enzyme to bind.
          </li>
          <li>
            <strong>BsmBI recognition site (pink):</strong> The sequence
            recognized by BsmBI enzyme.
          </li>
          <li>
            <strong>6nt sequence with 4bp overhang (blue):</strong> The sequence
            that will form the sticky ends after BsmBI cuts the DNA.
          </li>
          <li>
            <strong>Binding region (light gray):</strong> The region that binds
            to the target DNA sequence.
          </li>
        </ul>
        <p>The red arrows indicate where BsmBI cuts the DNA.</p>
      </div>
    </div>
  );
}

export default PrimerAnatomy;
