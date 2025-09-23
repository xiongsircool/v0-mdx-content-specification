import type { ReactNode } from "react"
import { Callout } from "@/components/mdx/callout"
import { Figure } from "@/components/mdx/figure"
import { CodeTabs } from "@/components/mdx/code-tabs"
import { Steps, Step } from "@/components/mdx/steps"
import { Badge } from "@/components/mdx/badge"

interface CalloutProps {
  type?: "info" | "success" | "warning" | "danger"
  children: ReactNode
}

// export function Callout({ type = "info", children }: CalloutProps) {
//   const variants = {
//     info: "border-blue-200 bg-blue-50 text-blue-900 dark:border-blue-800 dark:bg-blue-950 dark:text-blue-100",
//     success: "border-green-200 bg-green-50 text-green-900 dark:border-green-800 dark:bg-green-950 dark:text-green-100",
//     warning:
//       "border-yellow-200 bg-yellow-50 text-yellow-900 dark:border-yellow-800 dark:bg-yellow-950 dark:text-yellow-100",
//     danger: "border-red-200 bg-red-50 text-red-900 dark:border-red-800 dark:bg-red-950 dark:text-red-100",
//   }

//   return <div className={cn("my-6 rounded-lg border p-4", variants[type])}>{children}</div>
// }

interface FigureProps {
  src: string
  caption?: string
  alt?: string
}

// export function Figure({ src, caption, alt }: FigureProps) {
//   return (
//     <figure className="my-8">
//       <img src={src || "/placeholder.svg"} alt={alt || caption || ""} className="w-full rounded-lg border bg-muted" />
//       {caption && <figcaption className="mt-2 text-center text-sm text-muted-foreground">{caption}</figcaption>}
//     </figure>
//   )
// }

interface CodeTabsProps {
  tabs: string[]
  children: ReactNode
}

// export function CodeTabs({ tabs, children }: CodeTabsProps) {
//   return (
//     <div className="my-6">
//       <div className="flex border-b border-border">
//         {tabs.map((tab, index) => (
//           <button
//             key={tab}
//             className={cn(
//               "px-4 py-2 text-sm font-medium transition-colors",
//               index === 0 ? "border-b-2 border-primary text-primary" : "text-muted-foreground hover:text-foreground",
//             )}
//           >
//             {tab}
//           </button>
//         ))}
//       </div>
//       <div className="mt-4">{children}</div>
//     </div>
//   )
// }

interface StepsProps {
  children: ReactNode
}

// export function Steps({ children }: StepsProps) {
//   return <ol className="my-6 ml-6 list-none space-y-4 [counter-reset:step]">{children}</ol>
// }

export const mdxComponents = {
  Callout,
  Figure,
  CodeTabs,
  Steps,
  Step,
  Badge,
}
